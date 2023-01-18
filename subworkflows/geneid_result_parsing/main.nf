 process joinGffs {
    // show in the log which input file is analysed
    tag "joining gffs to ${output_file}"

    input:
    path (gffs_outputs)
    val (output_file)

    output:
    path ("${output_file}")

    script:
    """
    cat ${gffs_outputs} | egrep -v '^# ' | egrep -vw '###' | sort -u | sort -k1,1 -k4,5n > ${output_file}
    rm ${gffs_outputs}
    """
}


/*
 * Adding information from matches
 */
process intersectHints {
    // indicates to use as a container the value indicated in the parameter
    container "quay.io/biocontainers/bedtools:2.27.1--he513fc3_4"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${annotations_file}"

    input:
    path (annotations_file)
    path (matches_file)

    output:
    path ("${annotations_file_name}.labelled.tsv")

    script:
    annotations_file_name = annotations_file.getSimpleName()
    """
    sort -k1,1 -k4,5n ${annotations_file} > ${annotations_file}.sorted;
    sort -k1,1 -k4,5n ${matches_file} > ${matches_file}.sorted;

    bedtools intersect -sorted -a ${annotations_file}.sorted -b ${matches_file}.sorted -loj > ${annotations_file_name}.intersected.tmp

    cut -f-9,18 ${annotations_file_name}.intersected.tmp >> ${annotations_file_name}.labelled.tsv
    """
}



process processLabels {

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    // show in the log which input file is analysed
    tag "${annotations_file}"

    input:
    path (annotations_file)

    output:
    path ("${annotations_file_name}.gff3")

    script:
    annotations_file_name = annotations_file.BaseName
    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    data = pd.read_table("${annotations_file}",
                         names=["seq", "source", "feature", "start", "end", "strand",
                               "score", "phase", "attributes", "protein"],
                         header = None)

    grouped_data = data.groupby( ["seq", "source", "feature", "start", "end",
                                  "strand", "score", "phase",
                                  "attributes"]).agg({"protein" : ','.join }).reset_index()

    updated_attributes = [ x if y == '.' else x + f";ProteinIDs={y}" for x, y in zip(grouped_data["attributes"].values, grouped_data["protein"].values) ]

    grouped_data["attributes"] = updated_attributes
    grouped_data.drop("protein", axis = 1, inplace = True)

    grouped_data.sort_values(["seq", "start", "end"],
                  ascending = [True, True, False], inplace = True)

    grouped_data.to_csv("${annotations_file_name}.gff3",
                           sep = "\t",
                           header = None,
                           index = None)
    """
}
/*
 *
 */
process splitGff3 {

    tag "${annotations_file}"

    input:
    path (annotations_file)

    output:
    path ("${annotations_file_name}.head")
    path ("${annotations_file_name}.content")

    script:
    annotations_file_name = annotations_file.getSimpleName()
    """
    egrep '^##' ${annotations_file} | sort -u > ${annotations_file_name}.head;
    egrep -v '^#' ${annotations_file} | sort -u | sort -k1,1 -k4,5n > ${annotations_file_name}.content;
    """
}


/*
 *
 */
process mergeGff3 {

    publishDir(params.OUTPUT, mode: "copy", pattern : "*.gff3")

    // show in the log which input file is analysed
    tag "${annotations_file_name}"

    input:
    path (annotations_file_head)
    path (annotations_file_content)

    output:
    path ("${annotations_file_name}.gff3")

    script:
    annotations_file_name = annotations_file_head.getSimpleName()
    """
    cat ${annotations_file_head} ${annotations_file_content} > ${annotations_file_name}.gff3;
    """
}

/*
 * zip and index of gff3 for the portal
 */
process indexGff3 {

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

workflow add_labels {

    take:
    annotation_file
    hints_file

    main:

    (gff3_head, gff3_content) = splitGff3(annotation_file)

    labelled_content = intersectHints(gff3_content, hints_file) | processLabels

    final_gff3 = mergeGff3(gff3_head, labelled_content)

    emit:
    final_gff3
}


workflow geneid_result_parsing {
    
    take:
    gff_files
    hsp_found
    output_name

    main:


    merged_gff = joinGffs(gff_files, output_name)

    labeled_gff = add_labels(merged_gff, hsp_found)

    (gff3, gff3_gz, gff3_gz_tbi) = indexGff3(labeled_gff)

    emit:
    gff3
    gff3_gz
    gff3_gz_tbi
}