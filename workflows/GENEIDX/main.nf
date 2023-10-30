
output_dir = "${params.outdir}"
protein_db_dir = "${output_dir}/proteins"
species_dir = "${output_dir}/species"
subwork_folder = "${projectDir}/subworkflows"

//imported subworkflows
include { geneid_param_creation } from "${subwork_folder}/geneid_param_creation"

include { parameter_file_creation } from "${subwork_folder}/parameter_file_creation"

include { uniref_query } from "${subwork_folder}/uniref_query" addParams(OUTPUT: protein_db_dir)

include { diamond_blastx } from "${subwork_folder}/diamond_blastx" addParams(OUTPUT: protein_db_dir)

include { genomic_regions_estimation } from "${subwork_folder}/genomic_regions_estimation"

include { geneid_execution } from "${subwork_folder}/geneid_execution"

include { geneid_result_parsing } from "${subwork_folder}/geneid_result_parsing" addParams(OUTPUT: species_dir)

workflow GENEIDX {
  
  take:
  genomes
  meta

  main:
  //returns tuple containing id, acceptor path, donor path, start path, stop path and stdout containing a comma separated list of param values 
  geneid_parameters = geneid_param_creation(meta)

  //tuple val(id), path(proteins_fasta) 
  uniref_proteins = uniref_query(meta)

  //tuple val(id), path(diamond_out) // check if proteins have been downloaded already
  hsp_files = diamond_blastx(genomes, uniref_proteins)

  //tuple val(id), path(initial), path(trans)
  combined_matrices =  genomic_regions_estimation(genomes, hsp_files) 

  //tuple val(id), path(param_file)
  param_files = parameter_file_creation(geneid_parameters, combined_matrices)
  
  geneid_inputs = param_files.join(hsp_files)

  //tuple val(id), path(predictions)
  predictions = geneid_execution(genomes, geneid_inputs)
  .collectFile(){ item -> [ "${item[0]}-.gff3", item[1].text ]}
  .map { it -> tuple(it.baseName.tokenize('-')[0] ,it)}
  
  //tuple val(id), val(taxid), path(predictions), path(hsp_found)
  results = geneid_result_parsing(meta, predictions, hsp_files)

  emit:
  results
}

