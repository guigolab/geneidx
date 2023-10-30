process REPEATMODELER {
    tag "$id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::repeatmodeler=2.0.2a" : null)
    container 'quay.io/biocontainers/repeatmodeler:2.0.2a--pl5321h9ee0642_1' 
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.2a--pl5321h9ee0642_1':
    //     'quay.io/biocontainers/repeatmodeler:2.0.2a--pl5321h9ee0642_1' }"

    input:
    tuple val(id), path(fasta)

    output:
    tuple val(id), path("consensi.fa")

    script:

    """
    BuildDatabase -name genome_source -engine ncbi $fasta
    RepeatModeler -engine ncbi -pa ${task.cpus} -database genome_source
    cp RM_*/consensi.fa .
    """
}
