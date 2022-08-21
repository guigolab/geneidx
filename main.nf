#!/usr/bin/env nextflow

nextflow.enable.dsl=2


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

// this prints the input parameters
log.info """
GENEID+BLASTx - NextflowPipeline
=============================================
output			: ${params.output}
genome			: ${params.genome}
prot_file		: ${params.prot_file}
taxon			: ${params.taxid}
"""
// param_file		: ${params.param_f}

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
    log.info 'This is the Geneid+BLASTx test pipeline in Nextflow'
    log.info 'Please define the genome file, the protein file,\n\t\tthe parameter files and the output!\n'
    log.info '\n'
    exit 1
}

/*
 * Defining the output folders.
 */
OutputFolder = "${params.output}"

OutputFolderInternal = "${OutputFolder}/internal"
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

genoom = file(params.genome)
// proteins_file = file(params.prot_file)


/*
 * Defining the module / subworkflow path, and include the elements
 */
subwork_folder = "${projectDir}/subworkflows"

include { UncompressFASTA } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolder)

include { prot_down_workflow } from "${subwork_folder}/getProteins" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { build_protein_DB } from "${subwork_folder}/build_dmnd_db" addParams(OUTPUT: OutputFolderInternal,
  LABEL:'fourcpus')

include { alignGenome_Proteins } from "${subwork_folder}/runDMND_BLASTx" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'fourcpus')

include { matchAssessment } from "${subwork_folder}/getTrainingSeq" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { creatingParamFile } from "${subwork_folder}/modifyParamFile" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

include { geneid_WORKFLOW } from "${subwork_folder}/geneid" addParams( LABEL:'singlecpu' )

include { prep_concat } from "${subwork_folder}/prepare_concatenation" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

include { concatenate_Outputs_once } from "${subwork_folder}/geneid_concatenate" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

include { gff3addInfo } from "${subwork_folder}/addMatchInfo" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

// compress and index fastas to be stored and published to the cluster
include { compress_n_indexFASTA } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

// compress and index gff3s to be stored and published to the cluster
include { gff34portal } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolderSpeciesTaxid)


/*
 * MAIN workflow definition.
 */
workflow {

  // Uncompress FASTA file in here so that I can provide it
  //    uncompressed to the downstream modules
  uncompressed_genome = UncompressFASTA(genoom)


  // none of the returned objects is used by downsteam processes
  compress_n_indexFASTA(uncompressed_genome)


  // if proteins_file provided use proteins file
  // else, use taxon to download the proteins_file
  // both conditions are evaluated inside the execution of this workflow
  if (params.prot_file) {
    proteins_file = file(params.prot_file)
  } else {
    proteins_file = prot_down_workflow(params.taxid)
  }

  // Build protein database for DIAMOND
  protDB = build_protein_DB(proteins_file)


  // Run DIAMOND to find matches between genome and proteins
  hsp_found = alignGenome_Proteins(protDB, uncompressed_genome)


  // Automatic computation of the parameter file
  new_mats = matchAssessment(uncompressed_genome, hsp_found,
                                  params.match_score_min,
                                  params.match_ORF_min,
                                  params.intron_margin,
                                  params.min_intron_size,
                                  params.max_intron_size)


  new_param = creatingParamFile(params.no_score,
                                params.site_factor,
                                params.exon_factor,
                                params.hsp_factor,
                                params.exon_weight,

                                params.min_intron_size_geneid,
                                params.max_intron_size_geneid,

                                params.start_pwm,
                                params.acceptor_pwm,
                                params.donor_pwm,
                                params.stop_pwm,

                                new_mats.ini_comb,
                                new_mats.trans_comb
                                )


  // Run Geneid
  predictions = geneid_WORKFLOW(uncompressed_genome, new_param, hsp_found)

  // Prepare concatenation
  // This will initialize the final GFF3 file
  output_file = prep_concat(proteins_file, genoom)

  // Run concatenation of individual GFF3 files
  final_output = concatenate_Outputs_once(predictions.collect(), output_file)

  // Add information about the proteins to the final GFF3
  labelled_output = gff3addInfo(final_output, hsp_found)

  // fix gff3 file and compress it for the portal
  // gff34portal(final_output)
  gff34portal(labelled_output)

}


/*
 *  When complete print a message
 */
workflow.onComplete {
	println ( workflow.success ? "\nDone!\n" : "Oops ...\n" )
}
