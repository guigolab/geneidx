process FASTASPLITTER {
    tag "$id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::multiqc:1.12" : null)
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0':
    //     'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(id), path(genome)
    val(fsize) // size of fasta chunks to produce

    output:
    tuple val(id), path("*.part-*.*")

    script:
    """
       fasta-splitter.pl -part-sequence-size $fsize $genome
    """
}
