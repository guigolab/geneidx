// WorkflowGenomeannotator.initialise(params, log)

repeats_folder = "${params.output}/repeats"

subwork_folder = "${projectDir}/subworkflows"
modules_folder = "${projectDir}/modules"


include { repeat_download } from "${subwork_folder}/repeat_download" addParams(OUTPUT: repeats_folder,
 LABEL:'singlecpu')

include { REPEATMASKER } from "${subwork_folder}/REPEATMASKER"
include { REPEATMODELER } from "${modules_folder}/local/repeatmodeler"

workflow GENOMEANNOTATOR {


  take:
  assemblies //tuple id, assembly path
  metadata // tuple id, taxid

  main:

  ch_repeats = repeat_download(metadata, file(params.repeats_data_path))
  
  masked_assemblies = REPEATMASKER(assemblies, ch_repeats)

  emit:
  masked_assemblies
}

