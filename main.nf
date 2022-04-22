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
output				: ${params.output}
genome				: ${params.genome}
prot_file			: ${params.prot_file}
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
proteins_file = file(params.prot_file)


/*
 * Defining the module / subworkflow path, and include the elements
 */
subwork_folder = "${projectDir}/subworkflows/"

include { UncompressFASTA } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolder)


include { filter_Fasta_by_length } from "${subwork_folder}/filter_fasta" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { build_protein_DB } from "${subwork_folder}/build_dmnd_db" addParams(OUTPUT: OutputFolderInternal,
  LABEL:'fourcpus')

include { alignGenome_Proteins } from "${subwork_folder}/runDMND_BLASTx" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'fourcpus')

include { geneid_WORKFLOW } from "${subwork_folder}/geneid" addParams( LABEL:'singlecpu' )

include { concatenate_Outputs } from "${subwork_folder}/geneid_concatenate" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { matchAssessment } from "${subwork_folder}/getTrainingSeq" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { creatingParamFile } from "${subwork_folder}/modifyParamFile" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

// compress and index fastas to be stored and published to the cluster
include { compress_n_indexFASTA } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolderSpeciesTaxid)

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

  // Remove contigs that are too small (user chooses threshold)
  // filtered_genome = filter_Fasta_by_length(uncompressed_genome, params.min_seq_length)

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
  main_database_name = proteins_file.BaseName.toString().replaceAll("\\.fa", "")
  main_genome_name = genoom.BaseName.toString().replaceAll("\\.fa", "")
  // This is the name of the final GFF3 file
  out_filename = "${main_genome_name}.-.${main_database_name}.gff3"

  // Create the path to the file
  output_file = file(OutputFolder + "/" + out_filename)

  // Run concatenation of individual GFF3 files
  final_output = concatenate_Outputs(predictions, output_file)

  // fix gff3 file and compress it for the portal
  gff34portal(final_output.last())

}


/*
 *  When complete print a message
 */
workflow.onComplete {
	println ( workflow.success ? "\nDone!\n" : "Oops ...\n" )
}
