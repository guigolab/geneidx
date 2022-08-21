/*
 * Adding information from matches
 */
process gff3intersectHints {

    // where to store the results and in which way
    // publishDir(params.OUTPUT, mode : 'copy')

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
    annotations_file_name = annotations_file.BaseName
    """
    sort -k1,1 -k4,5n ${annotations_file} > ${annotations_file}.sorted;
    sort -k1,1 -k4,5n ${matches_file} > ${matches_file}.sorted;

    bedtools intersect -sorted -a ${annotations_file}.sorted -b ${matches_file}.sorted -loj > ${annotations_file_name}.intersected.tmp

    cut -f-9,18 ${annotations_file_name}.intersected.tmp >> ${annotations_file_name}.labelled.tsv
    """

}



process processLabels {

    // indicates to use as a container the value indicated in the parameter
    container "ferriolcalvet/python-modules"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${annotations_file}"

    input:
    path (annotations_file)

    output:
    path ("${annotations_file_name}.gff3")

    script:
    annotations_file_name = annotations_file.BaseName
    """
    #!/usr/bin/env python

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
process manageGff3sectionSplit {

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${annotations_file}"

    input:
    path (annotations_file)

    output:
    path ("${annotations_file_name}.head"), emit: head
    path ("${annotations_file_name}.content"), emit: content

    script:
    annotations_file_name = annotations_file.BaseName
    """
    egrep '^##' ${annotations_file} | sort -u > ${annotations_file_name}.head;
    egrep -v '^#' ${annotations_file} | sort -u | sort -k1,1 -k4,5n > ${annotations_file_name}.content;
    """
}


/*
 *
 */
process manageGff3sectionMerge {

    // where to store the results and in which way
    // publishDir(params.OUTPUT, mode : 'copy')

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${annotations_file_name}"

    input:
    path (annotations_file_head)
    path (annotations_file_content)

    output:
    path ("${annotations_file_name}.gff3")

    script:
    annotations_file_name = annotations_file_head.BaseName
    """
    cat ${annotations_file_head} ${annotations_file_content} > ${annotations_file_name}.gff3;
    """
}





/*
 * Workflow for adding information about the protein matches
 */

workflow gff3addInfo {

    take:
    annotation_file
    hints_file

    main:
    gff3_splitted = manageGff3sectionSplit(annotation_file)

    tsvWithLabels = gff3intersectHints(gff3_splitted.content, hints_file)
    gff3WithLabels = processLabels(tsvWithLabels)

    gff3_complete_labels = manageGff3sectionMerge(gff3_splitted.head, gff3WithLabels)

    emit:
    gff3_complete_labels

}
