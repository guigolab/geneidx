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

user_defined_params = params.grep( param -> param_keys.contains(param.key))


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
include { unzipFasta; createParamFile } from "${subwork_folder}/tools"
//imported subworkflows
include { geneid_param_creation } from "${subwork_folder}/geneid_param_creation"

include { diamond_blastx } from "${subwork_folder}/diamond_blastx" addParams(OUTPUT: protein_db_dir)

include { genomic_regions_estimation } from "${subwork_folder}/genomic_regions_estimation"

include { geneid_execution } from "${subwork_folder}/geneid_execution"

include { geneid_result_parsing } from "${subwork_folder}/geneid_result_parsing" addParams(OUTPUT: species_dir)

workflow {

  if(params.genome){
    genomes = channel.fromList([tuple("${params.taxid}_${params.result_name}", params.taxid, file(params.genome))])
  }else {
    genomes = channel.fromPath(params.tsv)
    .splitCsv( sep: '\t', header:true )
    .map { row -> tuple(row[params.row_id], row[params.row_taxid], file(row[params.row_path])) }
  }


  // send genomes from tsv


  // genome = channel.fromPath(params.genome)

  // if(!params.result_name){
  //   params.result_name = "${params.taxid}"
  // }
  // // create parameters to run geneid
  // (acc_pwm, don_pwm, sta_pwm, sto_pwm, geneid_param_values) = geneid_param_creation(params.taxid)

  // // Run DIAMOND to find matches between genome and proteins
  // hsp_found = diamond_blastx(genome)

  // // estimate genomic regions
  // new_mats = genomic_regions_estimation(genome, hsp_found)

  // // create parameter file to run geneid
  // new_param_file = createParamFile(
  //                                   params.taxid,
  //                                   geneid_param_values,
  //                                   sta_pwm,
  //                                   acc_pwm,
  //                                   don_pwm,
  //                                   sto_pwm,
  //                                   new_mats.ini_comb,
  //                                   new_mats.trans_comb,
  //                                   params.general_gene_params
  //                                 )

  // // execute geneid
  // predictions = geneid_execution(genome, new_param_file, hsp_found)

  // if (params.output_name) {
  //   output_name = param.result_file_name
  // } else {
  //   output_name = "${params.taxid}.gff3"
  // }

  //   (gff3, gff3_gz, gff3_gz_tbi) = geneid_result_parsing(predictions.collect(), hsp_found, output_name)

}
/*
 *  When complete print a message
 */
workflow.onComplete {
	println ( workflow.success ? "\nDone!\n" : "Oops ...\n" )
}
