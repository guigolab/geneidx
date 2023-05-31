subwork_folder = "${projectDir}/subworkflows/"

include { unzipFasta; getFASTA; getFASTA as getFASTA2} from "${subwork_folder}/tools"

/*
 * Remove some matches from the GFF to make it smaller and avoid redundancy in the introns
 */
process summarizeMatches {

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    // show in the log which input file is analysed
    tag "${name}"

    input:
    path (main_matches)

    output:
    path ("${name}.modified_exons.gff")

    script:
    name = main_matches.BaseName
    exon_margin = params.intron_margin
    """
    sort -k1,1 -k4,5n -k9,9 ${main_matches} | \
              awk '!found[\$1"\t"\$2"\t"\$3"\t"\$4]++' | \
              awk '!found[\$1"\t"\$2"\t"\$3"\t"\$5]++' > ${name}.summarized_matches

    awk 'OFS="\t"{print \$1, \$4-${exon_margin}, \$5+${exon_margin}, \$9\$7}' \
          ${name}.summarized_matches | \
          sort -k1,1 -k4,4 -k2,2n > ${name}.modified_exons.gff

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
    tag "${name}"

    input:
    path (main_matches)


    output:
    path ("${name}.introns.gff")

    script:
    min_intron_size=params.min_intron_size
    max_intron_size=params.max_intron_size
    name = main_matches.BaseName

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

    dd_intron.to_csv("${name}.introns.gff",
                     sep = "\t",
                     header = None,
                     index = None)
    """
}

/*
 * Use bedtools to remove the introns overlapping protein matches
 */
process removeProtOverlappingIntrons {

    // indicates to use as a container the value indicated in the parameter
    container "quay.io/biocontainers/bedtools:2.27.1--he513fc3_4"

    // indicates to use as a label the value indicated in the parameter
    // show in the log which input file is analysed
    tag "${introns_name}"

    input:
    path (main_matches)
    path (introns)

    output:
    path ("${introns_name}.non_overlapping_matches.gff")

    script:
    main_matches_name ="main_matches"
    introns_name = "introns"
    """
    sort -k1,1 -k4,5n ${introns} > ${introns}.sorted;
    sort -k1,1 -k4,5n ${main_matches} > ${main_matches}.sorted;

    bedtools intersect -sorted -a ${introns}.sorted \
                       -b ${main_matches}.sorted \
                       -v > ${introns_name}.non_overlapping_matches.gff
    """
}

/*
 * Merge the DIAMOND BLASTx matches and correct the scores
 * also change format to GFF3
 THIS COULD BE EASILY PARALLELIZABLE BUT I AM NOT SURE
 WHAT IS THE BEST WAY TO DO IT IN NEXTFLOW

 Requirement: blast2gff docker
 */
process mergeMatches {

    // show in the log which input file is analysed
    tag "${name}"

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    input:
    path(file_path)

    output:
    path ("${name}.gff3")

    script:
    name = file_path.BaseName
    
    """
    # get the sequences that have matches
    cut -f1 ${file_path} | uniq | sort -u > ${name}.seqs

    # iterate the sequences with matches, running blast2gff for each of them
    while read seq; do
        awk -F '\t' -v myvar=\$seq '\$1==myvar {print > (myvar".gff"); next} {print > ("REST.txt")}' ${file_path};
        mv REST.txt ${file_path};
        blast2gff -vg \${seq}.gff >> ${name}.SR.gff;
        rm \${seq}.gff;
    #    break
    done < ${name}.seqs;

    # remove the header rows of all files
    grep -v '#' ${name}.SR.gff > ${name}.SR.gff.tmp;
    mv ${name}.SR.gff.tmp ${name}.SR.gff;

    # change from GFF to GFF3 making each line a different "transcript"
    grep -v '#' ${name}.SR.gff | \
            awk '{\$3="CDS";print \$0,"ID="NR";Parent="NR";"}' | \
            sed -e 's/ /\t/g' > ${name}.gff3

    """
}

/*
 * Filter HSPs GFF3 by the score of the match
 */
process filterByScore {

    // show in the log which input file is analysed
    tag "${name}"

    // indicates to use as a label the value indicated in the parameter
    input:
    path(gff3_file)

    output:
    path ("${name}.over${score}.gff3")

    script:
    score = "${params.match_score_min}"
    name = gff3_file.BaseName
    """
    awk '\$6>=${score}' ${gff3_file} > ${name}.over${score}.gff3
    """
}

/*
 * Find ORFs
 */
process findORF {

    // indicates to use as a container the value indicated in the parameter
    container "quay.io/biocontainers/orfipy:0.0.4--py38h8c62d01_0"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${name}"

    input:
    path(seqs_file)

    output:
    path ("${name}_longest.bed")

    script:
    name = seqs_file.BaseName
    min_orf_len = params.match_ORF_min
    """
    orfipy ${seqs_file} --strand f --bed ${name}.bed --longest \
                        --min ${min_orf_len} --between-stops \
                        --outdir .
    """
}

/*
 * ORFs relative coordinates to absolute

 and have the python script in the container also
 see docker with "sgp_getHSPSR.pl" file for an example on how to do it.
 */
process updateCoords {

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    // show in the log which input file is analysed
    tag "${name}"

    input:
    path (original_gff3)
    path (relative_coords_file)

    output:
    path ("${name}.ORFs.gff3")

    script:
    name = original_gff3.BaseName
    """
    #!/usr/bin/env python3

    import pandas as pd

    read_orf_coords = pd.read_csv("${relative_coords_file}", sep='\t', header=None)
    read_gff3_file = pd.read_csv("${original_gff3}", sep='\t', header=None)

    read_orf_coords.columns = ["id", "rel_start", "rel_end", "info", "frame2", "strand_ORF"]
    read_gff3_file.columns = ["chr", "program", "region", "start", "end", "value", "strand", "frame", "id"]
    read_orf_coords.loc[:,"id"] = 'ID='+ read_orf_coords.loc[:,"id"].astype(str) + ';Parent='+ read_orf_coords.loc[:,"id"].astype(str) + ';'

    merged_orf_gff3 = read_gff3_file.merge(read_orf_coords, on='id', how='inner')

    def fix_coords(x):
        if x["strand"] == '+':
            return ( int(x["start"] + x["rel_start"]), int(x["start"] + x["rel_end"]) - 1 )

        return ( int(x["end"] - x ["rel_end"])+ 1, int(x["end"] - x["rel_start"]) )

    fixed_coords = merged_orf_gff3.apply(fix_coords, axis=1)
    fixed_coords = pd.DataFrame(fixed_coords.to_list())
    # print(fixed_coords)

    merged_orf_gff3.loc[:,"start"] = fixed_coords.iloc[:,0]
    merged_orf_gff3.loc[:,"end"]   = fixed_coords.iloc[:,1]

    # print(merged_orf_gff3.columns)
    ORF_coords = merged_orf_gff3[['chr', 'program', 'region', 'start', 'end', 'value', 'strand', 'frame', 'id']]
    ORF_coords.to_csv("${name}.ORFs.gff3", sep='\t', index=False, header=False)
    """
}

/*
 * Get the initial and transition probability matrices of the CDS
 */
process getCDSMatrices {

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    // show in the log which input file is analysed
    tag "${name}"

    input:
    path (cds)

    output:
    path(initial), emit: initial
    path(transition), emit: transition

    script:
    name = "${cds.BaseName}.cds"
    initial = "${name}.5.initial"
    transition = "${name}.5.transition"
    """
    FastaToTbl ${cds} > ${name}.tbl

    gawk '{print \$1,substr(\$2,1,length(\$2)-3)}' ${name}.tbl | \
                    gawk -f /scripts/MarkovMatrices.awk 5 ${name}

    sort +1 -2  -o ${initial} ${initial}
    sort +1 -2  -o ${transition} ${transition}
    """
}


/*
 * Get the initial and transition probability matrices of the introns
 */
process getIntronMatrices {

    // indicates to use as a container the value indicated in the parameter

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    // show in the log which input file is analysed
    tag "${introns_name}"

    input:
    path (introns)

    output:
    path (initial), emit: initial
    path (transition), emit: transition

    script:
    introns_name = "${introns.BaseName}.intron"
    initial = "${introns_name}.5.initial"
    transition = "${introns_name}.5.transition"
    """
    FastaToTbl ${introns} > ${introns_name}.tbl

    gawk -f /scripts/MarkovMatrices-noframe.awk 5 ${introns_name} ${introns_name}.tbl

    sort +1 -2  -o ${initial} ${initial}
    sort +1 -2  -o ${transition} ${transition}
    """
}


/*
 * Get the initial probability matrices of the introns
 */
process combineIni {

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    // show in the log which input file is analysed
    tag "${name}"

    input:
    path (cds_mats_ini)
    path (intron_mats_ini)

    output:
    path ("${name}.cds-intron.5.initial.geneid")

    script:
    name = "${cds_mats_ini.BaseName}_${intron_mats_ini.BaseName}"
    """
    ##  Compute log-likelihood exon matrices, assuming intron
    ##  matrices describing background probabilities

    gawk -f /scripts/pro2log_ini.awk ${intron_mats_ini} ${cds_mats_ini} \
          >  ${name}.cds-intron.5.initial

    gawk 'BEGIN {p=-1}{if (((NR+2) % 3)==0) p+=1; print \$2,p,\$1,\$3}' \
      ${name}.cds-intron.5.initial > ${name}.cds-intron.5.initial.geneid

    sed -i '1 i\\Markov_Initial_probability_matrix' ${name}.cds-intron.5.initial.geneid
    """
}

/*
 * Get the transition probability matrices of the introns
 */
process combineTrans {

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    // show in the log which input file is analysed
    tag "${name}"

    input:
    path (cds_mats_trans)
    path (intron_mats_trans)

    output:
    path ("${name}.cds-intron.5.transition.geneid")

    script:
    name = "${cds_mats_trans}_${intron_mats_trans}"
    """
    ##  Compute log-likelihood exon matrices, assuming intron
    ##  matrices describing background probabilities

    gawk -f /scripts/pro2log_tran.awk ${intron_mats_trans} ${cds_mats_trans} \
          >  ${name}.cds-intron.5.transition

    gawk 'BEGIN {p=-1}{if (((NR+2) % 3)==0) p+=1; print \$2,p,\$1,\$4}' \
      ${name}.cds-intron.5.transition > ${name}.cds-intron.5.transition.geneid

    sed -i '1 i\\Markov_Transition_probability_matrix' ${name}.cds-intron.5.transition.geneid
    """
}

/*
 * Workflow for obtaining the estimates of the exon sequences
 */
workflow cds_workflow {

    take:
    unzipped_genomes
    hsp_files

    main:

    hsp_gff3 = mergeMatches(hsp_files)

    filtered_hsp_gff3 = filterByScore(hsp_gff3)

    hsp_rel = getFASTA(filtered_hsp_gff3, unzipped_genomes) | findORF

    hsp_abs = updateCoords(hsp_gff3, hsp_rel)

    cds_matrix = getFASTA2(hsp_abs, unzipped_genomes) | getCDSMatrices

    emit:
    cds_mats_ini=cds_matrix.initial
    cds_mats_trans=cds_matrix.transition

}

/*
 * Workflow for obtaining the estimates of the intron sequences
 */

workflow intron_workflow {

    take:
    unzipped_genomes
    hsp_files

    main:

    computed_introns = summarizeMatches(hsp_files) | pyComputeIntrons


    non_overlapping_introns = removeProtOverlappingIntrons(hsp_files, computed_introns)


    intron_matrix = getFASTA(non_overlapping_introns, unzipped_genomes) | getIntronMatrices

    emit:
    intron_mats_ini=intron_matrix.initial
    intron_mats_trans=intron_matrix.transition

}


/*
 * Workflow connecting the different pieces
 */
workflow genomic_regions_estimation {

    // definition of input
    take:
    genomes
    hsp_files

    main:

    
    // unzipped_genomes = unzipFasta(genomes)

    // requirements:
    (cds_mats_ini, cds_mats_trans) = cds_workflow(genomes, hsp_files)

    (intron_mats_ini, intron_mats_trans) = intron_workflow(genomes, hsp_files)

    combine_matrices_ini = combineIni(cds_mats_ini, intron_mats_ini)
    combine_matrices_trans = combineTrans(cds_mats_trans, intron_mats_trans)


    emit:
    ini_comb = combine_matrices_ini
    trans_comb = combine_matrices_trans
}