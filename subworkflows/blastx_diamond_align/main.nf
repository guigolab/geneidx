
process runDIAMOND_getHSPs_GFF {

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy', pattern : '*.gff')

    label "diamond"

    // show in the log which input file is analysed
    tag "${main_genome_name} against ${main_database_name}"

    input:
    path(dmnd_database_file)
    path(reference_genome_file)

    output:
    path ("${main_genome_name}.${main_database_name}.hsp.gff")

    script:
    main_database_name = dmnd_database_file.BaseName
    main_genome_name = reference_genome_file.BaseName
    """
    if [ ! -s ${params.OUTPUT}/${main_genome_name}.${main_database_name}.hsp.gff ]; then
        echo "Running matches ${main_genome_name}.${main_database_name}"
        fmt6_custom='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qstrand qframe'
        diamond blastx --db ${dmnd_database_file} \
                       --query ${reference_genome_file} \
                       --max-target-seqs 0 \
                       --max-hsps 0  \
                       --outfmt \$fmt6_custom \
                       --evalue 0.0001 \
                       --block-size 0.8 \
                       --threads ${task.cpus} \
                       --out ${main_genome_name}.${main_database_name}.hsp.out

        awk 'BEGIN{OFS="\t"}{if (\$13=="-"){frame=-(\$14);print \$1,"blastx","hsp",\$8,\$7,\$12,\$13,frame,\$2}else if(\$13=="+"){print \$1,"blastx","hsp",\$7,\$8,\$12,\$13,\$14,\$2}}' \
                    ${main_genome_name}.${main_database_name}.hsp.out | sort -k1,1 -k4,5n -k9,9 \
                    | egrep -vw '^MT' | cut -d '|' -f1 > ${main_genome_name}.${main_database_name}.hsp.gff

    else
        echo "${main_genome_name}.${main_database_name}.hsp.gff already built"
        ln -s ${params.OUTPUT}/${main_genome_name}.${main_database_name}.hsp.gff ${main_genome_name}.${main_database_name}.hsp.gff;
    fi
    rm ${reference_genome_file}
    """
}




/*
 * Workflow connecting the different pieces
 */

workflow blastx_diamond_align {

    // definition of input
    take:
    prot_DB_file
    genome_file

    main:

    // directly generate the GFF file
    //  it also checks if the GFF file is already present in the output
    matchesGFF = runDIAMOND_getHSPs_GFF(prot_DB_file, genome_file)

    emit:
    matchesGFF

}
