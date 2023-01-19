/*
 * Uncompressing if needed
 */
process unzipFasta {

    label "geneidx"

    // show in the log which input file is analysed
    tag "${ref_to_index}"

    input:
    path (ref_to_index)

    output:
    path ("${main_genome_file}")

    script:
    main_genome_file = ref_to_index.BaseName
    print ref_to_index

    """
    if [ ! -s  ${main_genome_file} ]; then
        echo "unzipping genome ${main_genome_file}.gz"
        gunzip -c ${main_genome_file}.gz > ${main_genome_file};
    fi
    """
}



/*
 * Fixing chromosome names if needed
 */
process changeChromosomeName {

    label "geneidx"

    // show in the log which input file is analysed
    tag "${ref_to_index}"

    input:
    file (ref_to_index)

    output:
    path ("${main_genome_file}.clean.fa")

    script:
    main_genome_file = ref_to_index.BaseName
    """
    paste -d ' ' <(FastaToTbl ${main_genome_file}.fa | cut -d ' ' -f1 | cut -d '|' -f3) <(FastaToTbl ${main_genome_file}.fa | cut -d ' ' -f2) | FastaToTbl > ${main_genome_file}.clean.fa
    """
}

/*
 * Indexing if needed
 *    with the exonerate fastaindex
 */

process indexFasta {

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    // show in the log which input file is analysed
    tag "${main_genome_file}"

    input:
    path main_genome_file

    output:
    path ("${main_genome_file}.i")

    script:
    """
    if [ ! -s  ${main_genome_file}.i ]; then
        echo "indexing genome ${main_genome_file}"
        fastaindex -f ${main_genome_file} -i ${main_genome_file}.i
    fi
    """
}


/*
 * Indexing if needed
 *      using samtools faidx
 */
process faidxFasta {

    // indicates to use as a label the value indicated in the parameter
    label "samtools"

    // show in the log which input file is analysed
    tag "${main_genome_file}"

    input:
    path (main_genome_file)

    output:
    path ("${main_genome_file}.fai")

    script:
    """
    if [ ! -s  ${main_genome_file}.fai ]; then
        echo "indexing genome ${main_genome_file}"
        samtools faidx -f ${main_genome_file}
    fi
    """
}



/*
 * Use gffread to get the sequence of the introns
 */
 process getFASTA {

     // indicates to use as a container the value indicated in the parameter
     container "quay.io/biocontainers/gffread:0.12.7--hd03093a_1"

     // show in the log which input file is analysed
     tag "${gff3_name}"

     // indicates to use as a label the value indicated in the parameter
     label (params.LABEL)

     input:
     tuple val(gff3_name), path (gff3_file)
     tuple val(genome_name),path(genome_path)

     output:
     path ("${gff3_name}.fa")

     script:
     """
     gffread -x ${gff3_name}.fa -g ${genome_path} ${gff3_file};
     """
}


process createParamFile {
    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)
    // // show in the log which input file is analysed
    // tag "${output_param}"
    input:

    val (taxid)
    val (params_list)
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
    para_vals = params_list.split(',').collectEntries { entry ->
                                                  def pair = entry.split(':')
                                                  [(pair.first()): pair.last()]
                                                }
    
    no_score = Double.parseDouble(para_vals.no_score)

    absolute_cutoff_exons = Double.parseDouble(para_vals.absolute_cutoff_exons)
    coding_cutoff_oligos = Double.parseDouble(para_vals.coding_cutoff_oligos)

    site_factor = Double.parseDouble(para_vals.site_factor)
    exon_factor = Double.parseDouble(para_vals.exon_factor)
    hsp_factor = Double.parseDouble(para_vals.hsp_factor)
    exon_weight = Double.parseDouble(para_vals.exon_weights)

    ini_exon_weight = exon_weight + 1
    ini_coding_cutoff_oligos = coding_cutoff_oligos + 5
    output_param = file(params.genome).BaseName.toString() + ".${params.match_score_min}.${params.match_ORF_min}.manually_created.param"

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
