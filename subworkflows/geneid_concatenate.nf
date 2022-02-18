process concatenate_Outputs {
/*

POTENTIAL PROBLEM, MUTLIPLE PROCESSES ACCESSING AND GENERATING THE SAME FILE

*/



    // where to store the results and in which way
    publishDir(params.OUTPUT, mode: "copy", pattern : "*.gff3")

    // indicates to use as container the value indicated in the parameter
    // container params.CONTAINER

    // show in the log which input file is analysed
    // tag "${ref}"
    tag "concatenate Geneid outputs"


    input:
    path (gffs_outputs)
    val species_name
    val assembly_name


    output:
    path ("${main_genome_file}.gff3")

    script:
    main_genome_file = "${species_name}.${assembly_name}"
    partial_filename = gffs_outputs.name
    // location_files = params.OUTPUT
    // cat ${gffs_outputs} | grep -v '#' | sort -k1,1 >> ${main_genome_file}.gff3
    """
    cat ${params.OUTPUT}/${main_genome_file}*.gff3 | grep -v '#' | sort -u | sort -k1,1 >> ${main_genome_file}.gff3
    rm ${params.OUTPUT}/${partial_filename}
    """
}
