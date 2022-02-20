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
ENA_project_name    : ${params.ENA_project_name}
output				: ${params.output}
species				: ${params.species}
genome				: ${params.genome}
output				: ${params.output}
"""

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
    log.info 'This is Geneid\'s test pipeline in Nextflow'
    log.info 'Please define the ENA PROJECT ID and output!\n'
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

include { list_files_to_download } from "${subwork_folder}/files_to_download" addParams(OUTPUT: OutputFolder)
include { parse_json } from "${subwork_folder}/files_to_download" addParams(OUTPUT: OutputFolder)


include { DownloadFASTA_fromID } from "${subwork_folder}/download_fasta" addParams(OUTPUT: OutputFolder)

// include { geneid_WORKFLOW_single } from "${subwork_folder}/geneid_single" addParams(OUTPUT: OutputFolder)

include { geneid_WORKFLOW } from "${subwork_folder}/geneid" addParams(OUTPUT: OutputFolder, LABEL:'twocpus')
// include { UncompressFASTA } from "${subwork_folder}/geneid" addParams(OUTPUT: geneid_OutputFolder, LABEL:'twocpus')
// include { Index } from "${subwork_folder}/geneid" addParams(OUTPUT: geneid_OutputFolder, LABEL:'twocpus')
// include { runGeneid_fetching } from "${subwork_folder}/geneid" addParams(OUTPUT: geneid_OutputFolder, LABEL:'twocpus')
// include { concatenate_Outputs } from "${subwork_folder}/geneid_concatenate" addParams(OUTPUT: geneid_OutputFolder, LABEL:'twocpus')

include { build_protein_DB } from "${subwork_folder}/build_dmnd_db" addParams(OUTPUT: OutputFolder, LABEL:'twocpus')
include { alignGenome_Proteins } from "${subwork_folder}/runDMND_BLASTx" addParams(OUTPUT: OutputFolder, LABEL:'twocpus')




/*
 * MAIN workflow definition.
 */
workflow {

  // channel.from(params.ENA_project_name) | list_files_to_download | parse_json | flatten | set {taxons_set}
  //
  // downloaded_genome = DownloadFASTA_fromID(taxons_set)
  protDB = build_protein_DB(proteins_file)

  hsp_found = alignGenome_Proteins(protDB, genoom)
  // we call the runGeneid_fetching module using the channel for the queries
  // predictions = geneid_WORKFLOW_single(genoom, paar)

}


/*
 *  When complete print a message
 */
workflow.onComplete {
	println ( workflow.success ? "\nDone!\n" : "Oops ...\n" )
}
