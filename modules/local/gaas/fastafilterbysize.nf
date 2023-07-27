process GAAS_FASTAFILTERBYSIZE {
    tag "$id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gaas=1.2.0" : null)
    container 'quay.io/biocontainers/gaas:1.2.0--pl526r35_0'
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/gaas:1.2.0--pl526r35_0':
    //     'quay.io/biocontainers/gaas:1.2.0--pl526r35_0' }"

    input:
    tuple val(id), path(fasta)
    val min_size

    output:
    tuple val(id), path(filtered)

    script:
    filtered = "${id}.filtered.fa"
    """
    gaas_fasta_filter_by_size.pl -f $fasta -s $min_size -o $filtered

    """
}
