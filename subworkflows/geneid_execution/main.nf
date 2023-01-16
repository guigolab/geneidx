
OutputFolder = "${params.output}"

subwork_folder = "${projectDir}/subworkflows/"

include { indexFasta } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolder)


process runGeneid {

    label "geneidx"

    tag "run Geneid ${query}"

    input:
    path(reference_genome_file)
    path(reference_genome_index)
    path(geneid_param)
    path(protein_matches)
    val query

    output:
    path ("${main_genome_file}.*.gff3")

    script:
    main_genome_file = reference_genome_file.getName()
    main_output_file = protein_matches.getName().toString().replaceAll(".hsp", "")
    query_curated = query
    """
    # prepare sequence
    fastafetch -f ${reference_genome_file} -i ${reference_genome_index} -q \"${query}\" > ${main_genome_file}.${query}

    # prepare evidence
    egrep -w \"^${query}\" ${protein_matches} > ${main_output_file}.${query}.hsp.gff

    # check if the evidence file is empty
    if [ ! -s ${main_output_file}.${query}.hsp.gff ];  then
        echo "No protein evidence";
    else
        echo "Great, there is protein evidence";
    fi

    blast2gff -vg ${main_output_file}.${query}.hsp.gff > ${main_output_file}.${query}.SR.gff
    sgp_getHSPSR.pl \"${query}\" < ${main_output_file}.${query}.SR.gff > ${main_output_file}.${query}.HSP_SR.gff

    rm ${main_output_file}.${query}.hsp.gff
    rm ${main_output_file}.${query}.SR.gff

    # run Geneid + protein evidence
    geneid -3P ${geneid_param} -S ${main_output_file}.${query}.HSP_SR.gff ${main_genome_file}.${query} \
                | sed -e 's/geneid_v1.4/geneidx/g' | egrep -v '^# ' | grep -vFw '###' \
                >> ${main_output_file}.${query}.gff3

    rm ${main_output_file}.${query}.HSP_SR.gff
    rm ${main_genome_file}.${query}
    """
}

/*
 * Workflow connecting the different pieces
 */

workflow geneid_execution {

    take:
    ref_file
    geneid_param
    hsp_file


    main:

    genome = channel.fromPath(ref_file)

    index_filename = indexFasta(genome)

    genome.splitFasta( record: [id: true] )
                   .map{ it.toString().tokenize(':]').get(1) }
                   .set{ ch }

    predictions = runGeneid(genome,
                                      index_filename,
                                      geneid_param,
                                      hsp_file,
                                      ch)


    emit:
    predictions
}
