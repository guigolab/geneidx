/*
 * Uncompressing if needed
 */
process UncompressFASTA {

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
    // perl -i -lane 'if (/^>/) { (\$id, \$chr)=\$_=~/^>([\\w|.]+)[\\s\\w]+, [\\w]+: (\\w+)/; print ">".\$chr} else {print}' ${main_genome_file}
    // perl -i -lane 'if (/^>/) { ($id, $chr)=$_=~/^>([\w|.]+)[\s\w]+, chromosome: (\w+)/; print ">".$chr} else {print}' ${main_genome_file}
}


/*
 * Indexing if needed
 *    with the exonerate fastaindex
 */

process Index_i {

    // indicates to use as a container the value indicated in the parameter
    container "ferriolcalvet/geneid-fetching"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

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

    // copy output file
    // publishDir

    // indicates to use as a container the value indicated in the parameter
    container "quay.io/biocontainers/samtools:1.15--h1170115_1"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

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
