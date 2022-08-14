process concatenate_Outputs_once {
    // where to store the results and in which way
    // publishDir(params.OUTPUT, mode: "copy", pattern : "*.gff3")

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "adding to ${output_file}"

    input:
    path (gffs_outputs)
    path (output_file)

    output:
    path ("${output_file}")

    script:
    """
    cat ${gffs_outputs} | egrep -v '^# ' | egrep -vw '###' | sort -u | sort -k1,1 -k4,5n > ${output_file}
    rm ${gffs_outputs}
    """
}
