#!/usr/bin/env nextflow

/*
 *   This file was adapted from a part of 'nextflow-io/elixir-workshop-21'
 *   A course whose materials were prepared by Luca Cozzuto <lucacozzuto@gmail.com>
 */


/*
 * This code enables the new Nextflow dsl (domain-specific language).
 */

nextflow.enable.dsl=2


/*
 * @authors
 * Ferriol Calvet <ferriol.calvet@crg.eu>
 * Emilio Righi <emilio.righi@crg.eu>
 */

/*
 * Input parameters: ENA project name and output
 * The configuration is in nextflow.config file
 * Params are stored in the params.config file
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
param_file		: ${params.param_f}
"""

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
    log.info 'This is the Geneid+BLASTx test pipeline in Nextflow'
    log.info 'Please define the genome file, the protein file,\n\t\tthe parameter file and the output!\n'
    log.info '\n'
    exit 1
}

/*
 * Defining the output folders.
 */
OutputFolder = "${params.output}"

genoom = file(params.genome)
paar = file(params.param_f)
proteins_file = file(params.prot_file)


/*
 * Defining the module / subworkflow path, and include the elements
 */
subwork_folder = "${projectDir}/subworkflows/"


include { build_protein_DB } from "${subwork_folder}/build_dmnd_db" addParams(OUTPUT: OutputFolder,
  LABEL:'fourcpus')

include { alignGenome_Proteins } from "${subwork_folder}/runDMND_BLASTx" addParams(OUTPUT: OutputFolder,
  LABEL:'fourcpus')

include { geneid_WORKFLOW } from "${subwork_folder}/geneid" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { concatenate_Outputs } from "${subwork_folder}/geneid_concatenate" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')



/*
 * MAIN workflow definition.
 */
workflow {

  // Build protein database for DIAMOND
  protDB = build_protein_DB(proteins_file)

  // Run DIAMOND to find matches between genome and proteins
  hsp_found = alignGenome_Proteins(protDB, genoom)

  // **TO DO**
  // Evaluate matches and ORFs in matches
  // First to get an idea of how good/bad we are getting those regions
  //    later on we will generate a parameter file and use that file
  //    for running Geneid

  // Run Geneid
  predictions = geneid_WORKFLOW(genoom, paar, hsp_found)

  // Prepare concatenation
  main_database_name = proteins_file.BaseName.toString().replaceAll(".fa", "")
  main_genome_name = genoom.BaseName.toString().replaceAll(".fa", "")

  // This is the name of the final GFF3 file
  out_filename = "${main_genome_name}.-.${main_database_name}.gff3"

  // Create the path to the file
  output_file = file(OutputFolder + "/" + out_filename)

  // Run concatenation of individual GFF3 files
  final_output = concatenate_Outputs(predictions, output_file)

  // add header

}


/*
 *  When complete print a message
 */
workflow.onComplete {
	println ( workflow.success ? "\nDone!\n" : "Oops ...\n" )
}
