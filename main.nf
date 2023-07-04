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

output_dir = "${params.output}"
protein_db_dir = "${output_dir}/proteins"
species_dir = "${output_dir}/species"
subwork_folder = "${projectDir}/subworkflows"

//imported processes
include { createParamFile } from "${subwork_folder}/tools"
//imported subworkflows
include { geneid_param_creation } from "${subwork_folder}/geneid_param_creation"

include { parameter_file_creation } from "${subwork_folder}/parameter_file_creation"

include { uniref_query } from "${subwork_folder}/uniref_query" addParams(OUTPUT: protein_db_dir)

include { diamond_blastx } from "${subwork_folder}/diamond_blastx" addParams(OUTPUT: protein_db_dir)

include { genomic_regions_estimation } from "${subwork_folder}/genomic_regions_estimation"

include { geneid_execution } from "${subwork_folder}/geneid_execution"

include { geneid_result_parsing } from "${subwork_folder}/geneid_result_parsing" addParams(OUTPUT: species_dir)

workflow {
  //convert genome(s) to factory channel of tuples
  // if(params.genome){
  //   genomes = channel.fromList([tuple("${params.taxid}_${params.result_name}", params.taxid, file(params.genome))])
  // }else {
    // genomes = channel.fromPath(params.tsv)
    // .splitCsv( sep: '\t', header:true )
    // .map { row -> tuple(row[params.row_id], row[params.row_taxid], file(row[params.row_path])) }
  // }
  // genomes = channel.fromList([
  //   tuple("${params.taxid}_${params.result_name}", params.taxid, file(params.genome)),
  //   ])

  tsv = channel.fromPath(params.tsv)
  
  /*
  NOTES:
    nextflow manages the ena url in a weird way, it seems to stream the file instead of downloading it
    Maybe we should use another approach like installing ncbi datasets cli and download the assemblies from there
  */
  //genomes channel
  genomes = tsv.splitCsv( sep: '\t', header:true )
  .collectFile(){ row -> ["${row[params.row_id]}-.fa.gz", file(row[params.row_path])]}
  .map { it -> tuple(it.baseName.tokenize('-')[0],it)}
  //metadata channel
  meta = tsv.splitCsv( sep: '\t', header:true )
  .map { row -> tuple(row[params.row_id], row[params.row_taxid]) }

  //returns tuple containing id, acceptor path, donor path, start path, stop path and stdout containing a comma separated list of param values 
  geneid_parameters = geneid_param_creation(meta)

  //tuple val(id), path(proteins_fasta)
  uniref_proteins = uniref_query(meta)

  //tuple val(id), path(diamond_out)
  hsp_files = diamond_blastx(genomes, uniref_proteins)

  //tuple val(id), path(initial), path(trans)
  combined_matrices =  genomic_regions_estimation(genomes, hsp_files) 

  //tuple val(id), path(param_file)
  param_files = parameter_file_creation(geneid_parameters, combined_matrices)
  
  geneid_inputs = param_files.join(hsp_files)

  //tuple val(id), path(predictions)
  predictions = geneid_execution(genomes, geneid_inputs)
  .collectFile(){ item -> [ "${item[0]}-.gff3", item[1].text ]}
  .map { it -> tuple(it.baseName.tokenize('-')[0] ,it)} | view
  
  //tuple val(id), val(taxid), path(predictions), path(hsp_found)
  geneid_result_parsing(meta, predictions, hsp_files)

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
