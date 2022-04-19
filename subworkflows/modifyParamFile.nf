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
    val (no_score)

    val (site_factor)
    val (exon_factor)
    val (hsp_factor)
    val (exon_weight)

    val (min_intron_size)
    val (max_intron_size)

    path (start_pwm)
    path (acceptor_pwm)
    path (donor_pwm)
    path (stop_pwm)

    path (initial_probability_matrix)
    path (transition_probability_matrix)

    output:
    path ("${output_param}")

    script:
    ini_exon_weight = exon_weight + 1
    output_param = file(params.genome).BaseName.toString() + ".${params.match_score_min}.${params.match_ORF_min}.manually_created.param"
    """
    cat <<EOF > ${output_param}
    # geneid parameter file: human, 1 isochores
    # NO SCORE: ${no_score}
    # SITE FACTOR: ${site_factor}
    # EXON FACTOR: ${exon_factor}
    # HSP FACTOR: ${hsp_factor}
    # EXON WEIGHT: ${exon_weight}
    #
    # MIN INTRON SIZE: ${min_intron_size}
    # MAX INTRON SIZE: ${max_intron_size}
    #
    # START PWM: ${start_pwm}
    # ACCEPTOR PWM: ${acceptor_pwm}
    # DONOR PWM: ${donor_pwm}
    # STOP PWM: ${stop_pwm}
    #
    # INITIAL PWM: ${initial_probability_matrix}
    # TRANSITION PWM: ${transition_probability_matrix}
    #
    # Comment lines must start with '#'

    # Non-homology -0.35
    NO_SCORE
    ${no_score}

    # Number of isochores
    number_of_isochores
    1

    # PARAMETERS FROM THE SINGLE ISOCHORE

    # %GC
    boundaries_of_isochore
    0  100

    # Exons score: cutoffs
    Absolute_cutoff_exons
    -15 -15 -15 -15

    Coding_cutoff_oligos
    -10 -15 -15 -15

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


    # GENE MODEL: Rules about gene assembling (GenAmic)
    General_Gene_Model
    # INTRAgenic connections
    First+:Internal+                Internal+:Terminal+             ${min_intron_size}:${max_intron_size} block
    Terminal-:Internal-             First-:Internal-                ${min_intron_size}:${max_intron_size} blockr
    First+                          Intron+                         1:1 block
    Internal+                       Intron+                         1:1 block
    Intron+                         Internal+                       1:1 block
    Intron+                         Terminal+                       1:1 block
    Intron-                         First-                          1:1 block
    Intron-                         Internal-                       1:1 block
    Internal-                       Intron-                         1:1 block
    Terminal-                       Intron-                         1:1 block
    # External features
    Promoter+                       First+:Single+                  50:4000
    Terminal+:Single+               aataaa+                         50:4000
    First-:Single-                  Promoter-                       50:4000
    aataaa-                         Single-:Terminal-               50:4000
    # INTERgenic conections
    aataaa+:Terminal+:Single+       Single+:First+:Promoter+        500:Infinity
    aataaa+:Terminal+:Single+       Single-:Terminal-:aataaa-       500:Infinity
    Promoter-:First-:Single-        Single+:First+:Promoter+        500:Infinity
    Promoter-:First-:Single-        Single-:Terminal-:aataaa-       500:Infinity
    EOF

    cat ${output_param}.end >> ${output_param}
    rm ${output_param}.end;
    """
}
