process REPEATMASKER_REPEATMASK {
    tag "$id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::repeatmasker=4.1.2.p1" : null)
    container 'quay.io/biocontainers/repeatmasker:4.1.2.p1--pl5321hdfd78af_1'

    input:
    tuple val(id), path(fasta), env(LIBDIR), path(rm_lib)

    output:
    tuple val(id), path("*.masked")

    script:

    base_name = fasta.getName()
    genome_rm = base_name + ".masked"
    rm_gff = base_name + ".out.gff"
    options = "-lib $rm_lib"
    """
    # RepeatMasker $options -gff -xsmall -q -nolow -pa ${task.cpus} $fasta  # this one if for softmasking
    RepeatMasker $options -gff -q -nolow -pa ${task.cpus} $fasta            # this one if for hardmasking
    test -f ${genome_rm} || cp $fasta $genome_rm && touch $rm_gff
    """
}
