process CAT_FASTA {

    tag "$id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container 'biocontainers/biocontainers:v1.2.0_cv1'

    input:
    tuple val(id), path(fastas)

    output:
    tuple val(id), path(merged_fasta)

    script:
    def prefix = "${id}." + task.ext.prefix ?: "${id}"
    merged_fasta = prefix + ".fa"

    """
    cat $fastas > $merged_fasta

    """
}
