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

// keys of the geneid params (path) that can be tuned
param_keys = ['Acceptor_profile', 'Donor_profile', 'Start_profile', 'Stop_profile']


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
species_dir = "${output_dir}/species/${params.taxid}"
subwork_folder = "${projectDir}/subworkflows"

//imported processes
include { createParamFile } from "${subwork_folder}/tools"
//imported subworkflows
include { geneid_param_creation } from "${subwork_folder}/geneid_param_creation"

include { diamond_blastx } from "${subwork_folder}/diamond_blastx" addParams(OUTPUT: protein_db_dir)

include { genomic_regions_estimation } from "${subwork_folder}/genomic_regions_estimation"

include { geneid_execution } from "${subwork_folder}/geneid_execution"

include { geneid_result_parsing } from "${subwork_folder}/geneid_result_parsing"

workflow {
  //convert genome(s) to factory channel of tuples
  // if(params.genome){
  //   genomes = channel.fromList([tuple("${params.taxid}_${params.result_name}", params.taxid, file(params.genome))])
  // }else {
    genomes = channel.fromPath(params.tsv)
    .splitCsv( sep: '\t', header:true )
    .map { row -> tuple(row[params.row_id], row[params.row_taxid], file(row[params.row_path])) }
  // }
  // genomes = channel.fromList([
  //   tuple("1983395_${params.result_name}", "1983395", file(params.genome)),
  //   tuple("254365_${params.result_name}", "254365", file(params.genome))
  //   ])

    // genomes = channel.fromPath(params.tsv)
    // .splitCsv( sep: '\t', header:true )
    // .map { row -> tuple(row[params.row_id], row[params.row_taxid], file(row[params.row_path])) }

  // create parameters to run geneid
  (acc_pwm,don_pwm,sta_pwm,sto_pwm, param_values) = geneid_param_creation(genomes)

  //run diamond blastx
  hsp_found = diamond_blastx(genomes)

  // estimate genomic regions
  new_mats =  genomic_regions_estimation(genomes, hsp_found) 

  // create parameter file to run geneid
  new_param_file = createParamFile(
                                    genomes,
                                    param_values,
                                    sta_pwm,
                                    acc_pwm,
                                    don_pwm,
                                    sto_pwm,
                                    new_mats.ini_comb,
                                    new_mats.trans_comb,
                                  )
  // execute geneid
  predictions = geneid_execution(genomes, new_param_file, hsp_found) | view

  // if (params.output_name) {
  //   output_name = param.result_file_name
  // } else {
  //   output_name = "${params.taxid}.gff3"
  // }

  //parse and merge geneid predictions
  // (gff3, gff3_gz, gff3_gz_tbi) = geneid_result_parsing(predictions, hsp_found)

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
