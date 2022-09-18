/*
*  Geneid module.
*/

// Parameter definitions
// params.CONTAINER = "ferriolcalvet/training-geneid-params"
// params.OUTPUT = "geneid_output"
// params.LABEL = ""

/*
 * Update an empty parameter file
 */
process creatingParamFile {
    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy', pattern : '*.param')

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // // show in the log which input file is analysed
    // tag "${output_param}"

    input:
    // fixed values
    // boundaries_of_isochore
    // maximum_number_of_donors_per_acceptor_site
    // number_of_isochores

    // NO_SCORE

    // Absolute_cutoff_exons
    // Coding_cutoff_oligos

    // Exon_factor
    // Site_factor
    // HSP_factor
    // Exon_weights

    // Acceptor_profile
    // Donor_profile
    // Start_profile
    // Stop_profile

    // Markov_Initial_probability_matrix
    // Markov_Transition_probability_matrix

    // General_Gene_Model
    val (taxid)

    val (no_score)

    val (absolute_cutoff_exons)
    val (coding_cutoff_oligos)

    val (site_factor)
    val (exon_factor)
    val (hsp_factor)
    val (exon_weight)

    path (start_pwm)
    path (acceptor_pwm)
    path (donor_pwm)
    path (stop_pwm)

    path (initial_probability_matrix)
    path (transition_probability_matrix)

    path (general_model_values)


    output:
    path ("${output_param}")

    script:
    ini_exon_weight = exon_weight + 1
    ini_coding_cutoff_oligos = coding_cutoff_oligos + 5
    output_param = file(params.genome).BaseName.toString() + ".${params.match_score_min}.${params.match_ORF_min}.manually_created.param"
    // # MIN INTRON SIZE: ${min_intron_size}
    // # MAX INTRON SIZE: ${max_intron_size}
    """
    cat <<EOF > ${output_param}
    # geneid parameter file for ${taxid}, 1 isochores
    # PARAMETERS FROM A SINGLE ISOCHORE
    # NO_SCORE: ${no_score}
    # Site_factor: ${site_factor}
    # Exon_factor: ${exon_factor}
    # HSP_factor: ${hsp_factor}
    # Exon_weight: ${exon_weight}
    # Absolute_cutoff_exons: ${absolute_cutoff_exons}
    # Coding_cutoff_oligos: ${coding_cutoff_oligos} # we add +5 to the first exon value
    #
    #
    # START PWM: ${start_pwm}
    # ACCEPTOR PWM: ${acceptor_pwm}
    # DONOR PWM: ${donor_pwm}
    # STOP PWM: ${stop_pwm}
    #
    # INITIAL PWM: ${initial_probability_matrix}
    # TRANSITION PWM: ${transition_probability_matrix}
    #
    # General Gene Model parameters: ${general_model_values}
    # Comment lines must start with '#'

    # Non-homology
    NO_SCORE
    ${no_score}

    # Number of isochores
    number_of_isochores
    1

    # %GC
    boundaries_of_isochore
    0  100

    # Exons score: cutoffs
    Absolute_cutoff_exons
    ${absolute_cutoff_exons} ${absolute_cutoff_exons} ${absolute_cutoff_exons} ${absolute_cutoff_exons}

    Coding_cutoff_oligos
    ${ini_coding_cutoff_oligos} ${coding_cutoff_oligos} ${coding_cutoff_oligos} ${coding_cutoff_oligos}


    # Exon score: factors
    Site_factor
    ${site_factor} ${site_factor} ${site_factor} ${site_factor}

    Exon_factor
    ${exon_factor} ${exon_factor} ${exon_factor} ${exon_factor}

    HSP_factor
    ${hsp_factor} ${hsp_factor} ${hsp_factor} ${hsp_factor}

    Exon_weights
    ${ini_exon_weight} ${exon_weight} ${exon_weight} ${exon_weight}


    # Site prediction: Position Weight Arrays
    # Lenght, Offset, Cutoff and order (Markov model)
    EOF


    cat ${start_pwm} >> ${output_param};
    cat ${acceptor_pwm} >> ${output_param};
    cat ${donor_pwm} >> ${output_param};
    cat ${stop_pwm} >> ${output_param};


    cat <<EOF > ${output_param}.middle
    # Exon prediction: Markov model
    # Initial probabilities at every codon position
    Markov_oligo_logs_file
    5
    EOF

    cat ${output_param}.middle >> ${output_param};
    rm ${output_param}.middle;

    cat ${initial_probability_matrix} >> ${output_param};

    echo "# Transition probabilities at every codon position" >> ${output_param};

    cat ${transition_probability_matrix} >> ${output_param};


    cat <<EOF > ${output_param}.end

    # Donors per acceptor to build exons
    maximum_number_of_donors_per_acceptor_site
    5

    EOF

    cat ${output_param}.end >> ${output_param};

    cat ${general_model_values} >> ${output_param};

    rm ${output_param}.end;
    """
}
