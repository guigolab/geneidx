
OutputFolder = "${params.output}"

subwork_folder = "${projectDir}/subworkflows/"


process runGeneid {

    label "geneidx"

    tag "run Geneid"

    input:
    path(seq_file_path)
    path(geneid_param)
    path(protein_matches)

    output:
    path ("${seq_file_name}.gff3")

    script:
    seq_file_name = "geneid_result"
    """
    # prepare evidence

    query=\$(awk 'sub(/^>/, "")' ${seq_file_path})

    egrep -w \"^\$query\" ${protein_matches} > ${seq_file_name}.hsp.gff

    # check if the evidence file is empty
    if [ ! -s ${seq_file_name}.hsp.gff ];  then
        echo "No protein evidence";
    else
        echo "Great, there is protein evidence";
    fi

    blast2gff -vg ${seq_file_name}.hsp.gff > ${seq_file_name}.SR.gff
    sgp_getHSPSR.pl \"${seq_file_name}\" < ${seq_file_name}.SR.gff > ${seq_file_name}.HSP_SR.gff

    # run Geneid + protein evidence
    geneid -3P ${geneid_param} -S ${seq_file_name}.HSP_SR.gff ${seq_file_path} \
                | sed -e 's/geneid_v1.4/geneidx/g' | egrep -v '^# ' | grep -vFw '###' \
                >> ${seq_file_name}.gff3

    rm ${seq_file_name}.HSP_SR.gff
    """
}

/*
 * Workflow connecting the different pieces
 */

workflow geneid_execution {

    take:
    genome
    geneid_param
    hsp_file


    main:

    sequence = genome.splitFasta( file: true )

    predictions = runGeneid(sequence, geneid_param, hsp_file)


    emit:
    predictions
}
