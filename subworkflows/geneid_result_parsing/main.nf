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
    tuple val(id), path(annotations_file), path(matches_file)

    output:
    tuple val(id),path("${annotations_file_name}.labelled.tsv")

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
    tuple val(id), path(annotations_file)

    output:
    tuple val(id), path("${annotations_file_name}.gff3")

    script:
    annotations_file_name = annotations_file.BaseName
    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    data = pd.read_table("${annotations_file}",
                         names=["seq", "source", "feature", "start", "end", "strand",
                               "score", "phase", "attributes", "protein"],
                         header = None,
                         dtype='unicode')

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

process splitGff3 {

    tag "${annotations_file}"

    input:
    tuple val(id), path(annotations_file)

    output:
    tuple val(id), path(head)
    tuple val(id), path(content)

    script:
    head = "${id}.head"
    content = "${id}.content"
    """
    egrep '^##' ${annotations_file} | sort -u > ${head};
    egrep -v '^#' ${annotations_file} | sort -u | sort -k1,1 -k4,5n > ${content};
    """
}

process mergeGff3 {

    // show in the log which input file is analysed
    tag "${annotations_file_name}"

    input:
    tuple val(id), path(annotations_file_head), path(annotations_file_content)

    output:
    tuple val(id), path("${annotations_file_name}.gff3")

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
    publishDir "${params.OUTPUT}/${taxid}", mode : 'copy'

    // indicates to use as a container the value indicated in the parameter
    label "samtools"

    // indicates to use as a label the value indicated in the parameter
    // label (params.LABEL)

    // show in the log which input file is analysed
    tag "${annotations_file}"

    input:
    tuple val(id), val(taxid), path(annotations_file)

    output:
    tuple val(id), path("${annotations_file}"), path("${annotations_file}.gz"), path("${annotations_file}.gz.tbi")

    script:
    """
    egrep '^##' ${annotations_file} | sort -u > ${annotations_file}.head;
    egrep -v '^#' ${annotations_file} | sort -u | sort -k1,1 -k4,5n > ${annotations_file}.content;

    cat ${annotations_file}.head ${annotations_file}.content <(echo '###') > ${annotations_file};

    echo "Compressing ${annotations_file}";
    bgzip -c ${annotations_file} > ${annotations_file}.gz;

    echo "Indexing ${annotations_file}.gz";
    tabix -p gff ${annotations_file}.gz;
    """
}

workflow geneid_result_parsing {
    
    take:
    meta
    predictions
    hsp_files

    main:

    (gff3_head, gff3_content) = splitGff3(predictions)

    labeled_content = intersectHints(gff3_content.join(hsp_files)) | processLabels

    labeled_predictions = mergeGff3(gff3_head.join(labeled_content))

    output = indexGff3(meta.join(labeled_predictions))

    emit:
    output
    // output
}