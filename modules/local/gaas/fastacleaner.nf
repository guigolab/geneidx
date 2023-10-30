process GAAS_FASTACLEANER {
    tag "$id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gaas=1.2.0" : null)
    container 'quay.io/biocontainers/gaas:1.2.0--pl526r35_0'

    input:
    tuple val(id), path(fasta)

    output:
    tuple val(id), path(fasta_clean)

    script:
    fasta_clean = "${id}.clean.fa"
    """
    sed "s/[.]\$//" $fasta | sed "s/ .*//" > tmp
    gaas_fasta_cleaner.pl -f tmp -o $fasta_clean
    rm tmp
    """
}
