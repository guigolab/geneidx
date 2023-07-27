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
nextflow.enable.dsl=2
params.help             = false

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

wk_folder = "${projectDir}/workflows"

include { GENEIDX } from "${wk_folder}/GENEIDX"

include { GENOMEANNOTATOR } from "${wk_folder}/GENOMEANNOTATOR"

workflow {

  //get params --> genomeannotator --> geneidx
  //convert genome(s) to factory channel of tuples
  if(params.genome){
    genomes = channel.fromList([tuple(params.taxid, params.taxid, file(params.genome))])
    metadata = channel.fromList([tuple(params.taxid, params.taxid)])
  
  }else {

    tsv = channel.fromPath(params.tsv)

    genomes = tsv.splitCsv( sep: '\t', header:true )
    .collectFile(){ row -> ["${row[params.row_id]}-.fa.gz", file(row[params.row_path])]}
    .map { it -> tuple(it.baseName.tokenize('-')[0], it)}

    metadata = tsv.splitCsv( sep: '\t', header:true )
    .map { row -> tuple(row[params.row_id], row[params.row_taxid]) }
    }

    masked_genomes = GENOMEANNOTATOR(genomes, metadata)


}
/*
 *  When complete print a message
 */
workflow.onComplete {
	println ( "\nDone!\n" )
}

// workflow.onError {
//   println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
// }
