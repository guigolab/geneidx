process UNZIP_FASTA {

    label "geneidx"
    // show in the log which input file is analysed
    tag "${id}"

    input:
    tuple val(id), path(genome)

    output:
    tuple val(id), path(unzipped_genome)

    script:
    unzipped_genome = "${id}_unzipped.fa"
    """
    gunzip -c ${genome} > ${unzipped_genome};
    """
}


include { GAAS_FASTACLEANER } from '../../modules/local/gaas/fastacleaner'
include { GAAS_FASTAFILTERBYSIZE } from '../../modules/local/gaas/fastafilterbysize'

workflow ASSEMBLY_PREPROCESS {
    take:
    assemblies // tuple val(id), path(genome)

    main:

    unzipped_assemblies = UNZIP_FASTA(assemblies)
       	
    filtered_assemblies = GAAS_FASTAFILTERBYSIZE(unzipped_assemblies, params.min_contig_size) | GAAS_FASTACLEANER

    emit:
    filtered_assemblies
}
