/*
*  Geneid module.
*/

// Parameter definitions
params.CONTAINER = "ferriolcalvet/training-geneid-params"
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

include { Index_fai } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolder)

include { intron_workflow } from "${subwork_folder}/introns_estimates" addParams(OUTPUT: OutputFolder)

include { cds_workflow } from "${subwork_folder}/CDS_estimates" addParams(OUTPUT: OutputFolder)




/*
 * Get the initial and transition probability matrices of the CDS
 */
process getCDS_matrices {

    // indicates to use as a container the value indicated in the parameter
    container params.CONTAINER

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${cds_name}"

    input:
    path (cds)

    output:
    path ("${cds_name}.5.initial"), emit: initial
    path ("${cds_name}.5.transition"), emit: transition

    script:
    cds_name = cds.BaseName
    """
    FastaToTbl ${cds} > ${cds_name}.tbl

    gawk '{print \$1,substr(\$2,1,length(\$2)-3)}' ${cds_name}.tbl | \
                    gawk -f /scripts/MarkovMatrices.awk 5 ${cds_name}

    sort +1 -2  -o ${cds_name}.5.initial ${cds_name}.5.initial
    sort +1 -2  -o ${cds_name}.5.transition ${cds_name}.5.transition

    rm ${cds_name}.tbl
    """
}





/*
 * Get the initial and transition probability matrices of the introns
 */
process getIntron_matrices {

    // indicates to use as a container the value indicated in the parameter
    container params.CONTAINER

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${introns_name}"

    input:
    path (introns)

    output:
    path ("${introns_name}.5.initial"), emit: initial
    path ("${introns_name}.5.transition"), emit: transition

    script:
    introns_name = introns.BaseName
    """
    FastaToTbl ${introns} > ${introns_name}.tbl

    gawk -f /scripts/MarkovMatrices-noframe.awk 5 ${introns_name} ${introns_name}.tbl

    sort +1 -2  -o ${introns_name}.5.initial ${introns_name}.5.initial
    sort +1 -2  -o ${introns_name}.5.transition ${introns_name}.5.transition

    rm ${introns_name}.tbl
    """
}




/*
 * Get the initial probability matrices of the introns
 */
process CombineIni {

    // indicates to use as a container the value indicated in the parameter
    container params.CONTAINER

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${cds_name}"

    input:
    path (cds_mats_ini)
    path (intron_mats_ini)

    output:
    path ("${cds_name}.cds-intron.5.initial.geneid")

    script:
    cds_name = cds_mats_ini.BaseName
    """
    ##  Compute log-likelihood exon matrices, assuming intron
    ##  matrices describing background probabilities

    gawk -f /scripts/pro2log_ini.awk ${intron_mats_ini} ${cds_mats_ini} \
          >  ${cds_name}.cds-intron.5.initial

    gawk 'BEGIN {p=-1}{if (((NR+2) % 3)==0) p+=1; print \$2,p,\$1,\$3}' \
      ${cds_name}.cds-intron.5.initial > ${cds_name}.cds-intron.5.initial.geneid

    sed -i '1 i\\Markov_Initial_probability_matrix' ${cds_name}.cds-intron.5.initial.geneid
    """
}




/*
 * Get the transition probability matrices of the introns
 */
process CombineTrans {

    // indicates to use as a container the value indicated in the parameter
    container params.CONTAINER

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${cds_name}"

    input:
    path (cds_mats_trans)
    path (intron_mats_trans)

    output:
    path ("${cds_name}.cds-intron.5.transition.geneid")

    script:
    cds_name = cds_mats_trans.BaseName
    """
    ##  Compute log-likelihood exon matrices, assuming intron
    ##  matrices describing background probabilities

    gawk -f /scripts/pro2log_tran.awk ${intron_mats_trans} ${cds_mats_trans} \
          >  ${cds_name}.cds-intron.5.transition

    gawk 'BEGIN {p=-1}{if (((NR+2) % 3)==0) p+=1; print \$2,p,\$1,\$4}' \
      ${cds_name}.cds-intron.5.transition > ${cds_name}.cds-intron.5.transition.geneid

    sed -i '1 i\\Markov_Transition_probability_matrix' ${cds_name}.cds-intron.5.transition.geneid
    """
}


/*
 * Workflow connecting the different pieces
 */
workflow matchAssessment {

    // definition of input
    take:
    ref_file
    geneid_param
    hsp_file
    min_match_score
    min_match_ORF
    intron_margins
    min_intron
    max_intron

    main:

    // requirements:
    ref_file_ind = Index_fai(ref_file)

    cds_seq = cds_workflow(ref_file, ref_file_ind, hsp_file,
                            min_match_score, min_match_ORF)
    cds_mats = getCDS_matrices(cds_seq)
    cds_mats_ini = cds_mats.initial
    cds_mats_trans = cds_mats.transition

    introns_seq = intron_workflow(ref_file, ref_file_ind, hsp_file,
                            intron_margins, min_intron, max_intron)
    intron_mats = getIntron_matrices(introns_seq)
    intron_mats_ini = intron_mats.initial
    intron_mats_trans = intron_mats.transition

    combine_matrices_ini = CombineIni(cds_mats_ini, intron_mats_ini)
    combine_matrices_trans = CombineTrans(cds_mats_trans, intron_mats_trans)


    emit:
    ini_comb = combine_matrices_ini
    trans_comb = combine_matrices_trans
}
