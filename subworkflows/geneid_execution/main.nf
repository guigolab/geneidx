process runGeneid {

    label "geneidx"

    tag "run Geneid against ${seq_file_name}"

    input:
    tuple val(id), path(seq_file_path), path(geneid_param), path(protein_matches)

    output:
    tuple val(id), path("${seq_file_name}.gff3")

    script:
    seq_file_name = seq_file_path.BaseName
    """
    # prepare evidence

    query=\$(awk 'sub(/^>/, "")' ${seq_file_path})

    egrep -w \"^\$query\" ${protein_matches} > ${seq_file_name}.hsp.gff

    blast2gff -vg ${seq_file_name}.hsp.gff > ${seq_file_name}.SR.gff
    sgp_getHSPSR.pl \"${seq_file_name}\" < ${seq_file_name}.SR.gff > ${seq_file_name}.HSP_SR.gff

    # run Geneid + protein evidence
    geneid -3P ${geneid_param} -S ${seq_file_name}.HSP_SR.gff ${seq_file_path} \
                | sed -e 's/geneid_v1.4/geneidx/g' | egrep -v '^# ' | grep -vFw '###' \
                >> ${seq_file_name}.gff3

    """
}

/*
 * Workflow connecting the different pieces
 */

workflow geneid_execution {

    take:
    genomes
    geneid_inputs


    main:

    args = genomes.splitFasta(file: true).combine(geneid_inputs, by: 0)

    predictions = runGeneid(args)

    emit:
    predictions
}
