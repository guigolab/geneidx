/*
 * Uncompressing if needed
 */
process UncompressFASTA {

    label "geneidx"

    // show in the log which input file is analysed
    tag "${ref_to_index}"

    input:
    file (ref_to_index)

    output:
    path ("${main_genome_file}")

    script:
    main_genome_file = ref_to_index.BaseName

    """
    if [ ! -s  ${main_genome_file} ]; then
        echo "unzipping genome ${main_genome_file}.gz"
        gunzip -c ${main_genome_file}.gz > ${main_genome_file};
    fi
    """
}



/*
 * Fixing chromosome names if needed
 */
process fix_chr_names {

    label "geneidx"

    // show in the log which input file is analysed
    tag "${ref_to_index}"

    input:
    file (ref_to_index)

    output:
    path ("${main_genome_file}.clean.fa")

    script:
    main_genome_file = ref_to_index.BaseName
    """
    paste -d ' ' <(FastaToTbl ${main_genome_file}.fa | cut -d ' ' -f1 | cut -d '|' -f3) <(FastaToTbl ${main_genome_file}.fa | cut -d ' ' -f2) | FastaToTbl > ${main_genome_file}.clean.fa
    """
}







/*
 * Indexing if needed
 *    with the exonerate fastaindex
 */

process Index_i {

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    // show in the log which input file is analysed
    tag "${main_genome_file}"

    input:
    path main_genome_file

    output:
    path ("${main_genome_file}.i")

    script:
    """
    if [ ! -s  ${main_genome_file}.i ]; then
        echo "indexing genome ${main_genome_file}"
        fastaindex -f ${main_genome_file} -i ${main_genome_file}.i
    fi
    """
}





/*
 * Indexing if needed
 *      using samtools faidx
 */
process Index_fai {

    // indicates to use as a container the value indicated in the parameter

    // indicates to use as a label the value indicated in the parameter
    label "samtools"

    // show in the log which input file is analysed
    tag "${main_genome_file}"

    input:
    path (main_genome_file)

    output:
    path ("${main_genome_file}.fai")

    script:
    """
    if [ ! -s  ${main_genome_file}.fai ]; then
        echo "indexing genome ${main_genome_file}"
        samtools faidx -f ${main_genome_file}
    fi
    """
}






/*
 * Compressing fasta for the portal
 */
process compress_n_indexFASTA {

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy')

    // indicates to use as a container the value indicated in the parameter
    label "samtools"

    // show in the log which input file is analysed
    tag "${genome_file}"

    input:
    file (genome_file)

    output:
    path ("${genome_file}.gz")
    path ("${genome_file}.gz.gzi")
    path ("${genome_file}.gz.fai")

    script:
    """
    echo "compressing ${genome_file}";
    bgzip -i ${genome_file};

    echo "Indexing ${genome_file}";
    samtools faidx ${genome_file}.gz;
    """
}

/*
 * Indexing for the portal
 */
process gff34portal {

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy')

    // indicates to use as a container the value indicated in the parameter
    label "samtools"

    // indicates to use as a label the value indicated in the parameter
    // label (params.LABEL)

    // show in the log which input file is analysed
    tag "${annotations_file}"

    input:
    path (annotations_file)

    output:
    path ("${annotations_file}")
    path ("${annotations_file}.gz")
    path ("${annotations_file}.gz.tbi")

    script:
    """
    egrep '^##' ${annotations_file} | sort -u > ${annotations_file}.head;
    egrep -v '^#' ${annotations_file} | sort -u | sort -k1,1 -k4,5n > ${annotations_file}.content;

    cat ${annotations_file}.head ${annotations_file}.content <(echo '###') > ${annotations_file};

    rm ${annotations_file}.head;
    rm ${annotations_file}.content;

    echo "Compressing ${annotations_file}";
    bgzip -c ${annotations_file} > ${annotations_file}.gz;

    echo "Indexing ${annotations_file}.gz";
    tabix -p gff ${annotations_file}.gz;
    """

}








/*
 * Use gffread to get the sequence of the introns
 */
 process getFASTA {

     // indicates to use as a container the value indicated in the parameter
     container "quay.io/biocontainers/gffread:0.12.7--hd03093a_1"

     // show in the log which input file is analysed
     tag "${gff3_file}"

     // indicates to use as a label the value indicated in the parameter
     label (params.LABEL)

     input:
     path (gff3_file)
     path (ref_genome)
     path (ref_genome_index)

     output:
     path ("${main_gff3_file}.fa")

     script:
     main_gff3_file = gff3_file.BaseName
     """
     gffread -x ${main_gff3_file}.fa -g ${ref_genome} ${gff3_file};
     """
}
