/*
*  Geneid module.
*/

// Parameter definitions
params.CONTAINER = "ferriolcalvet/training-modules"
// params.OUTPUT = "geneid_output"
// params.LABEL = ""






/*
 * Merge the DIAMOND BLASTx matches and correct the scores
 * also change format to GFF3
 */
process mergeMatches {

    // // indicates to use as a container the value indicated in the parameter
    // build the docker container with blast2gff alone
    // container params.CONTAINER

    // show in the log which input file is analysed
    tag "${gff3_file}"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    input:
    file (gff3_file)

    output:
    path ("${main_genome_file}.fa")

    script:
    main_genome_file = gff3_file.BaseName

    """
    blast2gff -vg ${main_output_file}.hsp.gff > ${main_output_file}.SR.gff

    // change format from GFF to GFF3
    grep -v '#' ${main_output_file}.SR.gff | \
                awk '{\$3="CDS";print \$0,"ID="NR";Parent="NR";"}' | \
                sed -e 's/ /\t/g' > ${main_output_file}.hsp.gff3
    """
}


/*
 * Filter HSPs GFF3 by the score of the match
 */
process filter_by_score {

    // // indicates to use as a container the value indicated in the parameter
    // build the docker container with blast2gff alone
    // container params.CONTAINER

    // show in the log which input file is analysed
    tag "${gff3_file}"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    input:
    file (gff3_file)
    val score

    output:
    path ("${main_gff3_file}.over${score}.gff3")

    script:
    main_gff3_file = gff3_file.BaseName
    """
    awk '\$6>=${score}' ${main_gff3_file} > ${main_gff3_file}.over${score}.gff3
    """
}




/*
 * Evaluate against a reference GFF3
 * Provide as output a TSV with the stats?
 */
process evaluateGFF3 {

    // indicates to use as a container the value indicated in the parameter
    container params.CONTAINER

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${main_genome_file}"

    input:
    path (main_genome_file)

    output:
    path ("${main_genome_file}.fai")

    script:
    """
    gffcompare -T -r ${reference_gff3} ${query_GFF3} -o ${query_GFF3}.out
    rm ${query_GFF3}.out.loci
    rm ${query_GFF3}.out.annotated.gtf
    rm ${query_GFF3}.out.combined.gtf
    rm ${query_GFF3}.out.tracking
    """
}


/*
 * Uncompressing if needed
 */
process UncompressFASTA {

    // show in the log which input file is analysed
    tag "${ref_to_index}"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    input:
    path (ref_to_index)

    output:
    path ("${main_genome_file}")

    script:
    main_genome_file = ref_to_index.BaseName

    """
    if [ ! -s  ${main_genome_file} ]; then
        echo "unzipping genome ${main_genome_file}.gz"
        gunzip -c ${main_genome_file}.gz > ${main_genome_file};
    fi
    """
}

/*
 * Get sequence of transcripts from GFF3 file
 */
process getFASTA {

    // // indicates to use as a container the value indicated in the parameter
    // container params.CONTAINER

    // show in the log which input file is analysed
    tag "${gff3_file}"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    input:
    path (gff3_file)
    path (ref_genome)
    path (ref_genome_index)

    output:
    path ("${main_gff3_file}.fa")

    script:
    main_gff3_file = gff3_file.BaseName
    """
    gffread -x ${main_gff3_file}.fa -g ${ref_genome} ${gff3_file};
    """
}



/*
 * Indexing if needed
 */
process Index {

    // indicates to use as a container the value indicated in the parameter
    container params.CONTAINER

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

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
 * Find ORFs
 */
process ORF_finder {

    // indicates to use as a container the value indicated in the parameter
    container params.CONTAINER

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${main_genome_file}"

    input:
    path (main_genome_file)
    val min_orf_len

    output:
    path ("${main_genome_file}.fai")

    script:
    """
    orfipy ${seqs_file} --strand f --bed ${seqs_file}.bed --longest \
                        --min ${min_orf_len} --between-stops
                        --outdir .
    rm ${seqs_file}.bed
    # remove log also
    ls -l *.log | head -1 | xargs -n 1 rm
    """
}

/*
 * ORFs relative coordinates to absolute
 */
process updateGFFcoords {

    // indicates to use as a container the value indicated in the parameter
    container params.CONTAINER

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${main_genome_file}"

    input:
    path (original_gff3)
    path (relative_coords_file)

    output:
    path ("${main_genome_file}.fai")

    script:
    original_gff3_basename = original_gff3.BaseName
    """
    # change relative to absolute coordinates
    python python_script ${relative_coords_file} ${original_gff3} ${original_gff3_basename}.ORFs.gff3
    """
}



/*
 * Compute the Hexamer frequencies of the coding regions
 */
process trainHexFreq {

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // indicates to use as container the value indicated in the parameter
    container params.CONTAINER

    // show in the log which input file is analysed
    // tag "${ref}"
    tag "run Geneid ${query}"

    input:
    path(seq_file)

    output:
    path ("${main_genome_file}.*.gff3")

    script:
    main_genome_file = reference_genome_file.BaseName
    main_output_file = protein_matches.BaseName.toString().replaceAll(".hsp", "")
    """
    #!/bin/bash

    ##########################
    # Building coding statistic (Markov model order 5)
    #Generating the initial and transition matrices
    # without stop codon
    #################
    #markov5
    ##################

    gawk '{print \$1,substr(\$2,1,length(\$2)-3)}' all.cds_filter1.tbl | gawk -f MarkovMatrices.awk 5 set1.cds

    sort +1 -2  -o set1.cds.5.initial set1.cds.5.initial
    sort +1 -2  -o set1.cds.5.transition set1.cds.5.transition

    gawk -f MarkovMatrices-noframe.awk 5  set1.intron all.intron_filter1.tbl

    sort -o set1.intron.5.initial set1.intron.5.initial
    sort -o  set1.intron.5.transition set1.intron.5.transition


    ##  Compute log-likelihood exon matrices, assuming intron
    ##  matrices describing background probabilities

    gawk -f pro2log_ini.awk set1.intron.5.initial set1.cds.5.initial \
      >  set1.cds-intron.5.initial
    gawk -f pro2log_tran.awk set1.intron.5.transition set1.cds.5.transition \
      >  set1.cds-intron.5.transition

    gawk 'BEGIN {p=-1}{if (((NR+2) % 3)==0) p+=1; print \$2,p,\$1,\$3}' \
      set1.cds-intron.5.initial > set1.cds-intron.5.initial.geneid
    gawk 'BEGIN {p=-1}{if (((NR+2) % 3)==0) p+=1; print \$2,p,\$1,\$4}' \
      set1.cds-intron.5.transition > set1.cds-intron.5.transition.geneid

    """
    // https://genome.crg.es/software/geneid/training.html
}



/*
 * Update the parameter file for geneid
 */
process updateParamFile {

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // indicates to use as container the value indicated in the parameter
    container params.CONTAINER

    // show in the log which input file is analysed
    // tag "${ref}"
    tag "run Geneid ${query}"

    input:
    path(hex_freq_file)
    path(geneid_param)

    output:
    path ("${main_genome_file}.*.gff3")

    script:
    main_genome_file = reference_genome_file.BaseName
    main_output_file = protein_matches.BaseName.toString().replaceAll(".hsp", "")
    query_curated = query
    // we used this before when we were not cleaning the fasta identifiers
    // query_curated = query.toString().tokenize('|').get(1)
    """
    blast2gff -vg ${main_output_file}.${query}.hsp.gff > ${main_output_file}.${query}.SR.gff
    """
    // $projectDir/scripts/sgp_getHSPSR.pl \"${query}\" < ${main_genome_file}.${query}.SR.gff > ${main_genome_file}.${query}.HSP_SR.gff
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


    main:

    // // requirements:
    // gffcompare only for version 0.11 I am using version 0.12.6 in the cluster, maybe create a new container
    // gffread quay.io/biocontainers/gffread:0.12.7--hd03093a_1
    // python + modules pandas and some others
    // orfipy  quay.io/biocontainers/orfipy:0.0.4--py38h4a32c8e_1
    // (samtools)
    // dependencies of the Geneid training part


    // get matches from previous steps

    // pass matches through blast2gff
    // change format from gff to GFF3
    original_HSP_gff3 = mergeMatches(hsp_file)

    // filter the HSPs by the score
    score = 200
    filtered_HSPs_gff3 = filter_by_score(original_HSP_gff3, score)

    // report first metrics using gffcompare
    metrics1 = evaluateGFF3(filtered_HSPs_gff3, reference_gff3)

    // use gffread for obtaining the sequences of the matches
    // if the genome has an index next to it gffread goes much faster
    // samtools faidx genome.fa
    reference_genome_unc = UncompressFASTA(reference_genome)
    reference_genome_ind = Index(reference_genome_unc)

    matches_seqs = getFASTA(filtered_HSPs, reference_genome_unc)


    // use orfipy for obtaining the longest ORFs of each match
    // use python script for computing the absolute coordinates
    hspORFs_coords = ORF_finder(matches_seqs)


    // evaluate the newly generated GFF3
    metrics2 = evaluateGFF3(hspORFs_coords, reference_gff3)


    // get the sequences of the ORFs using gffread again
    hspORFs_seqs = getFASTA(filtered_HSPs, reference_genome_unc)



    // compute the hexamer frequencies and update the values in the parameter file
    // trainHexFreq
    // updateParamFile




    genome_filename = UncompressFASTA(ref_file)
    // genome_filename.subscribe {  println "Got: $it"  }

    index_filename = Index(genome_filename)
    // index_filename.subscribe {  println "Got: $it"  }

    genome_filename.splitFasta( record: [id: true] )
                   // .subscribe {  println "Got: $it"  }
                   .map{x -> x.toString().tokenize(':]').get(1)}
                   .set{ch}
                   // .subscribe {  println "Got: $it"  }
                   // .flatten()

    // ch.view()
    // ch.subscribe {  println "Got: $it"  }

    // genome_filename.view()
    // index_filename.view()

    // we call the runGeneid_fetching module using the channel for the queries
    predictions = runGeneid_fetching(genome_filename,
                                      index_filename,
                                      geneid_param,
                                      hsp_file,
                                      ch)


    emit:
    predictions
}
