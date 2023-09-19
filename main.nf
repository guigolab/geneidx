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

nextflow.enable.dsl=2

log.info """
GENOMEMASKING+GENEID+BLASTx - NextflowPipeline
=============================================
output path: ${params.output}
tsv input path: ${params.tsv}
id column: ${params.column_id_value}
taxid column: ${params.column_taxid_value}
path column: ${params.column_path_value}
mask genomes? ${params.use_masking}
"""

wk_folder = "${projectDir}/workflows"   
subwk_folder = "${projectDir}/subworkflows"

include { GENEIDX } from "${wk_folder}/GENEIDX"
include { GENOMEANNOTATOR } from "${wk_folder}/GENOMEANNOTATOR"
include { UNZIP_FASTA } from "${subwk_folder}/ASSEMBLY_PREPROCESS"

workflow {

    tsv = channel.fromPath(params.tsv)
    
    input_genomes = tsv.splitCsv( sep: '\t', header:true )
    .collectFile(){ row -> ["${row[params.column_id_value]}-.fa.gz", file(row[params.column_path_value])]}
    .map { it -> tuple(it.baseName.tokenize('-')[0], it)} | UNZIP_FASTA

    metadata = tsv.splitCsv( sep: '\t', header:true )
    .map { row -> tuple(row[params.column_id_value], row[params.column_taxid_value]) }

    genomes = input_genomes

    if(params.use_masking) genomes = GENOMEANNOTATOR(input_genomes, metadata)

    results = GENEIDX(genomes, metadata)

}


workflow.onComplete {
	println ( "\nDone!\n" )
}
