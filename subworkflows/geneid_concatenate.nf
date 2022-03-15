process concatenate_Outputs {
    /*
    POTENTIAL PROBLEM, MUTLIPLE PROCESSES ACCESSING AND GENERATING THE SAME FILE
    */

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode: "copy", pattern : "*.gff3")

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // indicates to use as container the value indicated in the parameter
    // container params.CONTAINER

    // show in the log which input file is analysed
    // tag "${ref}"
    tag "concatenate Geneid outputs"


    input:
    path (gffs_outputs)
    path (output_file)

    output:
    path ("${output_file}")

    script:
    // main_genome_file = matches_file.BaseName.toString().replaceAll(".hsp", "")
    // partial_filename = gffs_outputs.name
    // cat ${gffs_outputs} | grep -v '#' | sort -k1,1 >> ${main_genome_file}.gff3
    """
    cat ${gffs_outputs} | grep -v '#' | sort -u | sort -k1,1 -k4,5n >> ${output_file}
    rm ${gffs_outputs}
    """
    // cat ${params.OUTPUT}/${main_genome_file}*.gff3 | grep -v '#' | sort -u | sort -k1,1 >> ${main_genome_file}.gff3
    // rm ${params.OUTPUT}/${partial_filename}
    // rm ${gffs_outputs}

}
