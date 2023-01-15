#!/usr/bin/env nextflow

/*
 * @authors
 * Ferriol Calvet <ferriol.calvet@crg.eu>
 * Emilio Righi <emilio.righi@crg.eu>
 */

/*
 * Input parameters: genome, protein evidences, parameter file,
 * additional values for the generation of the parameter file.
 * Params are stored in the params.config file
 *
 * The configuration is in nextflow.config file
 */

// this prevents a warning of undefined parameter

params.help             = false

// keys of the geneid params (path) that can be tuned
param_keys = ['Acceptor_profile', 'Donor_profile', 'Start_profile', 'Stop_profile']

user_defined_params = params.grep( param -> param_keys.contains(param.key))
// this prints the input parameters
log.info """
GENEID+BLASTx - NextflowPipeline
=============================================
output			: ${params.output}
genome			: ${params.genome}
taxon			: ${params.taxid}
"""
// param_file		: ${params.param_f}

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
    log.info 'This is the GeneidX test pipeline in Nextflow'
    log.info 'Please define the genome file and the taxid of the species.\n'
    log.info 'To define additional parameters checkout the params.config file.\n'
    log.info '\n'
    exit 1
}

/*
 * Defining the output folders.
 */
OutputFolder = "${params.output}"

OutputFolderProteinDBs = "${OutputFolder}/proteins"
OutputFolderSpecies = "${OutputFolder}/species"
OutputFolderSpeciesTaxid = "${OutputFolder}/species/${params.taxid}"
// fasta.gz
// fasta.gz.gzi
// fasta.gz.fai
// gff3.gz
// gff3.gz.tbi
// hsp.gff
// param

// paramOutputFolder = "${params.output}/params"


/*
 * Defining the module / subworkflow path, and include the elements
 */
subwork_folder = "${projectDir}/subworkflows"

include { alignGenome_Proteins } from "${subwork_folder}/runDMND_BLASTx" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'fourcpus')

include { geneid_param_selection } from "${subwork_folder}/geneid_param_selection" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { geneid_param_profiles } from "${subwork_folder}/geneid_param_profiles" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { geneid_param_values } from "${subwork_folder}/geneid_param_values" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { uniprot_fasta_download } from "${subwork_folder}/uniprot_fasta_download" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu') 

include { diamond_db_build } from "${subwork_folder}/diamond_db_build" addParams(OUTPUT: OutputFolder,
  LABEL:'fourcpus') 

include { blastx_diamond_align } from "${subwork_folder}/blastx_diamond_align" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'fourcpus')

include { genomic_regions_estimation } from "${subwork_folder}/genomic_regions_estimation" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { creatingParamFile_frommap; creatingParamFile } from "${subwork_folder}/modifyParamFile" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

include { geneid_WORKFLOW } from "${subwork_folder}/geneid" addParams( LABEL:'singlecpu' )

include { prep_concat } from "${subwork_folder}/prepare_concatenation" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

include { concatenate_Outputs_once } from "${subwork_folder}/geneid_concatenate" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

include { gff3addInfo } from "${subwork_folder}/addMatchInfo" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

// compress and index fastas to be stored and published to the cluster
include { compress_n_indexFASTA; UncompressFASTA;  gff34portal } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

// compress and index gff3s to be stored and published to the cluster

param_file_input = params.geneid_param_file
/*
 * MAIN workflow definition.
 */
workflow {

    //split genome into sequence objects (id and sequence) 
    //nextflow manages the genome source (ex: http/file) and the unzip
  // genome = file(params.genome)
  // genome_path =  channel.fromPath(params.genome)

  genome = file(params.genome)

  param_file = param_file_input ? channel.fromPath(param_file_input) : geneid_param_selection(params.taxid).param_file

  (acc_pwm, don_pwm, sta_pwm, sto_pwm) = geneid_param_profiles(param_file)
  
  //TODO: check if all param profiles have been defined otherwise generate the missing onse
//   if (params.acceptor_pwm) 
// {    (acc_pwm, don_pwm, sta_pwm, sto_pwm) = params.acceptor_pwm, params.donor_pwm, params.start_pwm, params.stop_pwm
// }  else 
// {    (acc_pwm, don_pwm, sta_pwm, sto_pwm) = geneid_param_parse(param_file)
// }
  geneid_param_values = geneid_param_values(param_file, params.maps_param_values)
  
  proteins_file = params.prot_file ? file(params.prot_file) : uniprot_fasta_download(params.taxid,params.proteins_lower_lim,params.proteins_upper_lim)

  // Build protein database for DIAMOND
  protDB = diamond_db_build(proteins_file)

  // Run DIAMOND to find matches between genome and proteins
  hsp_found = blastx_diamond_align(protDB, genome)

  unzipped_genome = UncompressFASTA(genome)

  new_mats = genomic_regions_estimation(unzipped_genome, hsp_found,
                                  params.match_score_min,
                                  params.match_ORF_min,
                                  params.intron_margin,
                                  params.min_intron_size,
                                  params.max_intron_size)

    // if sites matrices provided, use them
  // else, use taxon to get the closest geneid param file
  
  new_param = creatingParamFile_frommap(
                                params.taxid,
                                geneid_param_values.params_map,
                                sta_pwm,
                                acc_pwm,
                                don_pwm,
                                sto_pwm,
                                new_mats.ini_comb,
                                new_mats.trans_comb,
                                params.general_gene_params
                                )

  // Run Geneid
  predictions = geneid_WORKFLOW(genome, new_param, hsp_found)

  // Prepare concatenation
  // This will initialize the final GFF3 file
  // output_file = prep_concat(proteins_file, genoom)

  // // Run concatenation of individual GFF3 files
  // final_output = concatenate_Outputs_once(predictions.collect(), output_file)

  // // Add information about the proteins to the final GFF3
  // labelled_output = gff3addInfo(final_output, hsp_found)

  // // If proteins from UniRef, add GO terms to the GFF3
  // if (params.source_uniprot)
  //   gff34portal(labelled_output)
  // else 
  //   gff34portal(labelled_output)
  

    // get proteins file by default value if not defined in params
  //   proteins_file = params.prot_file ? 
  //     file(params.prot_file) 
  //     : proteins_file = prot_down_workflow(params.taxid,  params.proteins_lower_lim,  params.proteins_upper_lim)

  // // Build protein database for DIAMOND
  //   protDB = build_protein_DB(proteins_file)

  // // Run DIAMOND to find matches between genome and proteins
  //   hsp_found = alignGenome_Proteins(protDB, sequences)

  // // Automatic computation of the parameter file
  //   new_mats = matchAssessment(sequences, hsp_found,
  //                                 params.match_score_min,
  //                                 params.match_ORF_min,
  //                                 params.intron_margin,
  //                                 params.min_intron_size,
  //                                 params.max_intron_size)

  // // get param file from closest taxon if not defined in params
  //   def param_keys = ['acceptor_pwm', 'donor_pwm', 'start_pwm', 'stop_pwm']

  // // check if all the values are defined in params
  //   def collected_params = params.grep( param -> param_keys.contains(param.key))

  // if( collected_params.size() == 0){
  //   param_file_sel = param_selection_workflow(params.taxid, 0, parameter_location)
  //   println param_file_sel
  // }
  // else if(collected_params.size() < param_keys.size()) {
  //   println 'Error'
  // } 
  // else {
  //   mapped_params = [*:collected_params]
  //   // param_file_sel = param_selection_workflow(params.taxid, 0, parameter_location)
  //   // acc_pwm = param_file_sel.acceptor_pwm
  //   // don_pwm = param_file_sel.donor_pwm
  //   // sta_pwm = param_file_sel.start_pwm
  //   // sto_pwm = param_file_sel.stop_pwm
  // }

  // para_vals = param_value_selection_workflow(params.taxid, 0,
  //                                           parameter_location,
  //                                           params.maps_param_values)

  // new_param = creatingParamFile_frommap(
  //                               params.taxid,

  //                               para_vals.params_map,

  //                               sta_pwm,
  //                               acc_pwm,
  //                               don_pwm,
  //                               sto_pwm,

  //                               new_mats.ini_comb,
  //                               new_mats.trans_comb,

  //                               params.general_gene_params

  //                               )

  // Run Geneid
  // predictions = geneid_WORKFLOW(uncompressed_sequences, new_param, hsp_found)

  // // Prepare concatenation
  // // This will initialize the final GFF3 file
  // output_file = prep_concat(proteins_file, params.genome)

  // // Run concatenation of individual GFF3 files
  // final_output = concatenate_Outputs_once(predictions.collect(), output_file)

  // // Add information about the proteins to the final GFF3
  // labelled_output = gff3addInfo(final_output, hsp_found)

  // // If proteins from UniRef, add GO terms to the GFF3
  // if (params.source_uniprot){
  //   // func_labelled_output = addGOterms(labelled_output)
  //   // gff34portal(func_labelled_output)
  //   gff34portal(labelled_output)
  // } else {
  //   // fix gff3 file and compress it for the portal
  //   gff34portal(labelled_output)
  // }

}


/*
 *  When complete print a message
 */
workflow.onComplete {
	println ( workflow.success ? "\nDone!\n" : "Oops ...\n" )
}
