/*
*  Geneid module.
*/

// Parameter definitions
params.CONTAINER = "ferriolcalvet/training-modules"
// params.OUTPUT = "geneid_output"
// params.LABEL = ""



/*
 * Defining the output folders.
 */
OutputFolder = "${params.output}"


/*
 * Defining the module / subworkflow path, and include the elements
 */
subwork_folder = "${projectDir}/subworkflows/"

include { UncompressFASTA } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolder)
include { Index_fai } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolder)
include { getFASTA } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolder)


/*
 * Remove some matches from the GFF to make it smaller and avoid redundancy in the introns
 */
process summarizeMatches {

    // // indicates to use as a container the value indicated in the parameter
    // container params.CONTAINER

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${main_matches_name}"

    input:
    path (main_matches)
    val exon_margin

    output:
    path ("${main_matches_name}.modified_exons.gff")

    script:
    main_matches_name = main_matches.BaseName
    """
    sort -k1,1 -k4,5n -k9,9 ${main_matches} | \
              awk '!found[\$1"\t"\$2"\t"\$3"\t"\$4]++' | \
              awk '!found[\$1"\t"\$2"\t"\$3"\t"\$5]++' > ${main_matches}.summarized_matches

    awk 'OFS="\t"{print \$1, \$4-${exon_margin}, \$5+${exon_margin}, \$9\$7}' \
          ${main_matches}.summarized_matches | \
          sort -k1,1 -k4,4 -k2,2n > ${main_matches_name}.modified_exons.gff

    rm ${main_matches}.summarized_matches
    """
    // awk -v exmar="${exon_margin}" 'OFS="\t"{print \$1, \$4-exmar, \$5+exmar, \$9\$7}' \
    //             ${main_matches}.summarized_matches | \
    //             sort -k1,1 -k4,4 -k2,2n > ${main_matches_name}.modified_exons.gff
}




/*
 * Use a python script for identifying the introns
 */
process pyComputeIntrons {

    // indicates to use as a container the value indicated in the parameter
    container "ferriolcalvet/python-modules"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${main_matches_name}"

    input:
    path (main_matches)
    val min_intron_size
    val max_intron_size

    output:
    path ("${main_matches_name}.introns.gff")

    script:
    main_matches_name = main_matches.BaseName

    """
    #!/usr/bin/env python

    import pandas as pd

    data = pd.read_csv("${main_matches}",
                       sep = "\t", header = None)

    prev_id = ""
    prev_end = 0
    intron_l = []

    for index, row in data.iterrows():
        if row[3] == prev_id:
            intron_l.append(list(row) + [prev_end])
        else:
            prev_id = row[3]
        prev_end = row[2]


    dd_intron = pd.DataFrame(intron_l)
    dd_intron = dd_intron[[0,4,1,3]]

    # separate the strand from the ID
    dd_intron[3] = dd_intron[3].apply(lambda x : x[-1])

    # filter out the introns that do not accomplish the conditions of size
    dd_intron = dd_intron[(dd_intron[1] - dd_intron[4]) > int(${min_intron_size})]
    dd_intron = dd_intron[(dd_intron[1] - dd_intron[4]) < int(${max_intron_size})].reset_index(drop = True).reset_index()

    dd_intron["dot"] = '.'
    dd_intron["source"] = 'hsp'
    dd_intron["type"] = 'CDS'

    # we use index column for providing an ID to each intron
    dd_intron["index"] = "ID=" + dd_intron["index"].astype(str) + ";Parent=" + dd_intron["index"].astype(str) + ";"

    # select columns and rename
    dd_intron = dd_intron[[0, 'source', 'type', 4, 1, 'dot', 3, 'dot', 'index']]
    dd_intron.columns = ['seq', 'source', 'type', 'start', 'end', 'score', 'strand', 'frame', 'ID']

    # change column type to integer for a proper sorting
    dd_intron['start'] = dd_intron['start'].astype(int)
    dd_intron['end'] = dd_intron['end'].astype(int)

    # remove duplicated introns if any
    dd_intron.drop_duplicates(subset = ['seq', 'source', 'type', 'start', 'end'],
                                    keep = "first", inplace = True)

    dd_intron = dd_intron.sort_values(by = ['seq', 'start', 'end']).reset_index(drop = True)

    dd_intron.to_csv("${main_matches_name}.introns.gff",
                     sep = "\t",
                     header = None,
                     index = None)
    """
    // awk '!found[\$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5]++' ${main_matches_name}.introns.gff | \
    //              sort -k1,1 -k4,5n > ${main_matches_name}.introns.non_redundant.gff
}




/*
 * Use bedtools to remove the introns overlapping protein matches
 */
process removeProtOverlappingIntrons {

    // indicates to use as a container the value indicated in the parameter
    container "quay.io/biocontainers/bedtools:2.27.1--he513fc3_4"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${introns_name}"

    input:
    path (main_matches)
    path (introns)

    output:
    path ("${introns_name}.non_overlapping_matches.gff")

    script:
    main_matches_name = main_matches.BaseName
    introns_name = introns.BaseName
    """
    sort -k1,1 -k4,5n ${introns} > ${introns}.sorted;
    sort -k1,1 -k4,5n ${main_matches} > ${main_matches}.sorted;

    bedtools intersect -sorted -a ${introns} \
                       -b ${main_matches} \
                       -v > ${introns_name}.non_overlapping_matches.gff
    """
}


/*
 * Workflow for obtaining the estimates of the intron sequences
 */

workflow intron_workflow {

    // definition of input
    take:
    ref_file
    ref_file_ind
    hsp_file


    main:
    // // requirements:
    // gffread quay.io/biocontainers/gffread:0.12.7--hd03093a_1
    // python + modules pandas and some others
    exon_margins = 40
    gff_reduced = summarizeMatches(hsp_file, exon_margins)

    min_size = 0
    max_size = 10000
    computed_introns = pyComputeIntrons(gff_reduced, min_size, max_size)

    non_overlapping_introns = removeProtOverlappingIntrons(hsp_file, computed_introns)

    introns_seq = getFASTA(non_overlapping_introns, ref_file, ref_file_ind)

    emit:
    introns_seq

}
