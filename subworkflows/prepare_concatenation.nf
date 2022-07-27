
process prep_concat {
    /*
    POTENTIAL PROBLEM, MUTLIPLE PROCESSES ACCESSING AND GENERATING THE SAME FILE
    Solution so far, if I generate the file in advance and I pass it as input
    to the process Nextflow manages it appropriately
    */

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode: "copy", pattern : "*.gff3")

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "create file for concatenation ${out_filename}"


    input:
    path (prot_file)
    path (geno_file)

    output:
    path ("${out_filename}")

    script:
    main_database_name = prot_file.BaseName.toString().replaceAll("\\.fa", "")
    // println(main_database_name)
    main_genome_name = geno_file.BaseName.toString().replaceAll("\\.fa", "")

    // This is the name of the final GFF3 file
    out_filename = "${main_genome_name}.-.${main_database_name}.gff3"
    // println(out_filename)

    // // Create the path to the file
    // output_file = file(params.OUTPUT + "/" + out_filename)
    """
    printf "" > ${out_filename};
    """
    // egrep -v '^# ' ${gffs_outputs} >> ${output_file}
    // cat ${gffs_outputs} | grep -v '#' | sort -u | sort -k1,1 -k4,5n >> ${output_file}
    // rm ${gffs_outputs}

}
