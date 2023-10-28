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
output path: ${params.outdir}
tsv input path: ${params.tsv}
id column: ${params.column_id_value}
taxid column: ${params.column_taxid_value}
path column: ${params.column_path_value}
mask genomes? ${params.use_masking}
"""

wk_folder = "${projectDir}/workflows"   
subwk_folder = "${projectDir}/subworkflows"
species_dir = "${params.outdir}/species"

include { GENEIDX } from "${wk_folder}/GENEIDX"
include { GENOMEANNOTATOR } from "${wk_folder}/GENOMEANNOTATOR"
include { ASSEMBLY_PREPROCESS } from "${subwk_folder}/ASSEMBLY_PREPROCESS"
include { UNZIP_FASTA } from "${subwk_folder}/ASSEMBLY_PREPROCESS"
include { parseFastaHeader } from "${subwk_folder}/tools"

workflow {

    tsv = channel.fromPath(params.tsv).splitCsv( sep: '\t', header:true )
    .collectFile(storeDir:params.assemblies_dir){ row -> ["${row[params.column_id_value]}-${row[params.column_taxid_value]}-.fa.gz", file(row[params.column_path_value])]}
    .map { it -> 
        def elements = it.baseName.tokenize('-')
        tuple(elements[0],elements[1], it)
         }
    .branch { 
        validFasta: it[2].countFasta() > 0
        invalidFasta: true 
        } 

    invalid_fasta = tsv.invalidFasta
    if(invalid_fasta.count()){
        invalid_fasta.count().view { it -> "A total of ${it} invalid assemblies have been found"}
        invalid_fasta.collectFile(storeDir:params.outdir){ item ->
            [ "INVALID_ASSEMBLIES.txt", "${item[0]}\t${item[1]}\t${item[2]}" + '\n' ]
        }
        .subscribe {
            println "Invalid assemblies have been saved in: $it"
        }
    } 

    assemblies = tsv.validFasta.map{it-> tuple(it[0], it[2])}

    metadata = tsv.validFasta.map { row -> tuple(row[0], row[1]) }

    filtered_assemblies = UNZIP_FASTA(assemblies) | ASSEMBLY_PREPROCESS | parseFastaHeader 

    target_genomes = filtered_assemblies

    if(params.use_masking) target_genomes = GENOMEANNOTATOR(filtered_assemblies, metadata)

    results = GENEIDX(target_genomes, metadata)

}


workflow.onComplete {
	println ( "\nDone!\n" )
}

