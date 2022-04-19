/*
*  Geneid module.
*/

// Parameter definitions
params.CONTAINER = "ferriolcalvet/geneid-fetching"
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
include { Index_i } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolder)


process runGeneid_fetching {

    // where to store the results and in which way
    // publishDir(params.OUTPUT, pattern : '*.gff3')

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // indicates to use as container the value indicated in the parameter
    container params.CONTAINER

    // show in the log which input file is analysed
    // tag "${ref}"
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
    main_genome_file = reference_genome_file.BaseName
    main_output_file = protein_matches.BaseName.toString().replaceAll(".hsp", "")
    query_curated = query
    // we used this before when we were not cleaning the fasta identifiers
    // query_curated = query.toString().tokenize('|').get(1)
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
                | sed -e 's/geneid_v1.4/geneidblastx/g' \
                >> ${main_output_file}.${query}.gff3


    rm ${main_output_file}.${query}.HSP_SR.gff

    rm ${main_genome_file}.${query}
    """
    // $projectDir/scripts/sgp_getHSPSR.pl \"${query}\" < ${main_genome_file}.${query}.SR.gff > ${main_genome_file}.${query}.HSP_SR.gff
    // geneid -3P ${geneid_param} -S ${main_output_file}.${query}.HSP_SR.gff ${main_genome_file}.${query} \
    //             | sed -e 's/geneid_v1.4/geneidblastx/g' | egrep 'CDS' | sort -k4,5n \
    //             >> ${main_output_file}.${query}.gff3
}


/*
 * Workflow connecting the different pieces
 */

workflow geneid_WORKFLOW {

    // definition of input
    take:
    ref_file
    geneid_param
    hsp_file


    main:

    index_filename = Index_i(ref_file)
    // index_filename.view()

    ref_file.splitFasta( record: [id: true] )
                   // .subscribe {  println "Got: $it"  }
                   .map{x -> x.toString().tokenize(':]').get(1)}
                   .set{ch}
                   // .subscribe {  println "Got: $it"  }
                   // .flatten()

    // we call the runGeneid_fetching module using the channel for the queries
    predictions = runGeneid_fetching(ref_file,
                                      index_filename,
                                      geneid_param,
                                      hsp_file,
                                      ch)


    emit:
    predictions
}
