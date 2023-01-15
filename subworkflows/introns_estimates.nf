/*
*
*/

// Parameter definitions
// params.CONTAINER = "ferriolcalvet/training-modules"
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
}




/*
 * Use a python script for identifying the introns
 */
process pyComputeIntrons {

    // indicates to use as a container the value indicated in the parameter
    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

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
    #!/usr/bin/env python3

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
}







/*
 * Use a python script for identifying the introns

process pyComputeIntrons_SS {

    // indicates to use as a container the value indicated in the parameter
    container "ferriolcalvet/geneidx"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${main_matches_name}"

    input:
    path (main_matches)
    val min_intron_size
    val max_intron_size

    output:
    path ("${main_matches_name}.introns.gff3")
    path ("${main_matches_name}.acceptors.gff3")
    path ("${main_matches_name}.donors.gff3")

    script:
    main_matches_name = main_matches.BaseName

    """
    #!/usr/bin/env python3

    import pandas as pd

    data = pd.read_csv("${main_matches}",
                       sep = "\t", header = None)

    prev_id = ""
    prev_end = 0
    intron_l = []


    def find_acceptor(x):
        # 24 intron and 3 exon
        if x.iloc[3] == "+":
            return (x.iloc[2] - 40, x.iloc[2] + 10)
        elif x.iloc[3] == "-":
            return (x.iloc[1] - 10, x.iloc[1] + 40)

    def find_donor(x):
        # 3 exon and 6 intron
        if x.iloc[3] == "+":
            return (x.iloc[1] - 10, x.iloc[1] + 20)
        elif x.iloc[3] == "-":
            return (x.iloc[2] - 20, x.iloc[2] + 10)


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
    dd_intron = dd_intron[(dd_intron[1] - dd_intron[4]) < int(${max_intron_size})].reset_index(drop = True)
    dd_intron.columns = ["seq", "start", "end", "strand"]

    # remove duplicated introns if any
    dd_intron.drop_duplicates(keep = "first", inplace = True)
    dd_intron.reset_index(drop = True, inplace = True)



    # compute acceptors and donors
    acceptor_df = pd.DataFrame(dd_intron.apply(find_acceptor, axis = 1).to_list())
    acceptor_df = pd.concat((acceptor_df, dd_intron), axis = 1)[["seq", 0, 1, "strand"]]
    acceptor_df.columns = ["seq", "start", "end", "strand"]

    donor_df = pd.DataFrame(dd_intron.apply(find_donor, axis = 1).to_list())
    donor_df = pd.concat((donor_df, dd_intron), axis = 1)[["seq", 0, 1, "strand"]]
    donor_df.columns = ["seq", "start", "end", "strand"]


    dd_intron["dot"] = '.'
    dd_intron["source"] = 'hsp_int'
    dd_intron["type"] = 'CDS'
    dd_intron = dd_intron.reset_index()

    acceptor_df["dot"] = '.'
    acceptor_df["source"] = 'hsp_acc'
    acceptor_df["type"] = 'CDS'
    acceptor_df = acceptor_df.reset_index()

    donor_df["dot"] = '.'
    donor_df["source"] = 'hsp_don'
    donor_df["type"] = 'CDS'
    donor_df = donor_df.reset_index()


    # we use index column for providing an ID to each intron, acceptor and donor
    # also select columns and rename
    dd_intron["index"] = "ID=" + dd_intron["index"].astype(str) + ";Parent=" + dd_intron["index"].astype(str) + ";"
    dd_intron = dd_intron[["seq", 'source', 'type', "start", "end", 'dot', "strand", 'dot', 'index']]

    acceptor_df["index"] = "ID=" + acceptor_df["index"].astype(str) + ";Parent=" + acceptor_df["index"].astype(str) + ";"
    acceptor_df = acceptor_df[["seq", 'source', 'type', "start", "end", 'dot', "strand", 'dot', 'index']]

    donor_df["index"] = "ID=" + donor_df["index"].astype(str) + ";Parent=" + donor_df["index"].astype(str) + ";"
    donor_df = donor_df[["seq", 'source', 'type', "start", "end", 'dot', "strand", 'dot', 'index']]


    # change column type to integer for a proper sorting
    dd_intron['start'] = dd_intron['start'].astype(int)
    dd_intron['end'] = dd_intron['end'].astype(int)


    dd_intron = dd_intron.sort_values(by = ['seq', 'start', 'end']).reset_index(drop = True)

    dd_intron.to_csv("${main_matches_name}.introns.gff3",
                     sep = "\t",
                     header = None,
                     index = None)

    acceptor_df.to_csv("${main_matches_name}.acceptors.gff3",
                     sep = "\t",
                     header = None,
                     index = None)

    donor_df.to_csv("${main_matches_name}.donors.gff3",
                     sep = "\t",
                     header = None,
                     index = None)
    """
    // awk '!found[\$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5]++' ${main_matches_name}.introns.gff | \
    //              sort -k1,1 -k4,5n > ${main_matches_name}.introns.non_redundant.gff
}

*/







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

    bedtools intersect -sorted -a ${introns}.sorted \
                       -b ${main_matches}.sorted \
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
    intron_margins
    min_size
    max_size


    main:
    // intron_margins = 40
    gff_reduced = summarizeMatches(hsp_file, intron_margins)

    // min_size = 0
    // max_size = 10000
    computed_introns = pyComputeIntrons(gff_reduced, min_size, max_size)

    non_overlapping_introns = removeProtOverlappingIntrons(hsp_file, computed_introns)

    introns_seq = getFASTA(non_overlapping_introns, ref_file, ref_file_ind)

    emit:
    introns_seq

}
