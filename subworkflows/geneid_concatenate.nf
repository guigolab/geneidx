
process concatenate_Outputs {
    /*
    POTENTIAL PROBLEM, MUTLIPLE PROCESSES ACCESSING AND GENERATING THE SAME FILE
    Solution so far, if I generate the file in advance and I pass it as input
    to the process Nextflow manages it appropriately
    */

    // where to store the results and in which way
    // publishDir(params.OUTPUT, mode: "copy", pattern : "*.gff3")

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "concatenating ${gffs_outputs}"


    input:
    path (gffs_outputs)
    path (output_file)

    output:
    path ("${output_file}")

    script:
    """
    egrep -v '^# ' ${gffs_outputs} | egrep -vw '###' | sort -u >> ${output_file}
    rm ${gffs_outputs}
    """
    // egrep -v '^# ' ${gffs_outputs} >> ${output_file}
    // cat ${gffs_outputs} | grep -v '#' | sort -u | sort -k1,1 -k4,5n >> ${output_file}
    // rm ${gffs_outputs}

}
