
process prep_concat {
    
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

    """
    printf "" > ${out_filename};
    """

}
