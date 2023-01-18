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

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
    log.info 'This is the GeneidX test pipeline in Nextflow'
    log.info 'Please define the genome file and the taxid of the species.\n'
    log.info 'To define additional parameters checkout the params.config file.\n'
    log.info '\n'
    exit 1
}

OutputFolder = "${params.output}"

OutputFolderProteinDBs = "${OutputFolder}/proteins"
OutputFolderSpeciesTaxid = "${OutputFolder}/species/${params.taxid}"
subwork_folder = "${projectDir}/subworkflows"

//imported processes
include { unzipFasta; createParamFile } from "${subwork_folder}/tools"
//imported subworkflows
include { geneid_param_creation } from "${subwork_folder}/geneid_param_creation"

include { uniprot_fasta_download } from "${subwork_folder}/uniprot_fasta_download" addParams(OUTPUT: OutputFolderProteinDBs)

include { diamond_db_build } from "${subwork_folder}/diamond_db_build" addParams(OUTPUT: OutputFolderProteinDBs)

include { blastx_diamond_align } from "${subwork_folder}/blastx_diamond_align" addParams(OUTPUT: OutputFolderSpeciesTaxid)

include { genomic_regions_estimation } from "${subwork_folder}/genomic_regions_estimation"

include { geneid_execution } from "${subwork_folder}/geneid_execution"

include { geneid_result_parsing } from "${subwork_folder}/geneid_result_parsing" addParams(OUTPUT: OutputFolderSpeciesTaxid)


param_file_input = params.geneid_param_file
/*
 * MAIN workflow definition.
 */

workflow {

  genome = file(params.genome)

  unzipped_genome = unzipFasta(genome) //unzip for genomic regions estimation workflow

  (acc_pwm, don_pwm, sta_pwm, sto_pwm, geneid_param_values) = geneid_param_creation(params.taxid)

  proteins_file = params.prot_file ? 
    file(params.prot_file) : 
    uniprot_fasta_download(params.taxid, params.proteins_lower_lim, 
                            params.proteins_upper_lim
                          )
  // Build protein database for DIAMOND
  prot_DB = diamond_db_build(proteins_file)

  // Run DIAMOND to find matches between genome and proteins
  hsp_found = blastx_diamond_align(prot_DB, unzipped_genome)

  new_mats = genomic_regions_estimation(unzipped_genome, hsp_found,
                                  params.match_score_min,
                                  params.match_ORF_min,
                                  params.intron_margin,
                                  params.min_intron_size,
                                  params.max_intron_size)

  new_param_file = createParamFile(
                                params.taxid,
                                geneid_param_values,
                                sta_pwm,
                                acc_pwm,
                                don_pwm,
                                sto_pwm,
                                new_mats.ini_comb,
                                new_mats.trans_comb,
                                params.general_gene_params
                              )

  predictions = geneid_execution(genome, new_param_file, hsp_found)

  output_name = "${params.taxid}_${genome.getSimpleName()}.gff3"

  (gff3, gff3_gz, gff3_gz_tbi) = geneid_result_parsing(predictions.collect(), hsp_found, output_name)

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

}
/*
 *  When complete print a message
 */
workflow.onComplete {
	println ( workflow.success ? "\nDone!\n" : "Oops ...\n" )
}
