// WorkflowGenomeannotator.initialise(params, log)

OutputFolder = "${params.output}"

subwork_folder = "${projectDir}/subworkflows"
modules_folder = "${projectDir}/modules"


include { repeat_download } from "${subwork_folder}/repeat_download" addParams(OUTPUT: OutputFolder,
 LABEL:'singlecpu')

include { ASSEMBLY_PREPROCESS } from "${subwork_folder}/ASSEMBLY_PREPROCESS"
include { REPEATMASKER } from "${subwork_folder}/REPEATMASKER"
include { REPEATMODELER } from "${modules_folder}/local/repeatmodeler"

workflow GENOMEANNOTATOR {


  take:
  assemblies //tuple id, assembly path
  metadata // tuple id, taxid

  main:

  filtered_assemblies = ASSEMBLY_PREPROCESS(assemblies)

  //tuple val(id), path(repeats)
  // ch_repeats = params.repeat_modeler ?  
  //   repeat_download(metadata, file(params.repeats_data_path)) :
  //   REPEATMODELER(filtered_assemblies) --> memory exceed? do we need to use this?


  ch_repeats = repeat_download(metadata, file(params.repeats_data_path))
  
  masked_assemblies = REPEATMASKER(filtered_assemblies, ch_repeats)
    
  //   ch_versions = ch_versions.mix(REPEATMASKER.out.versions)
  //   ch_genome_rm = REPEATMASKER.out.fasta

  emit:
  masked_assemblies
}


workflow.onComplete {
    NfcoreTemplate.summary(workflow, params, log)
}
