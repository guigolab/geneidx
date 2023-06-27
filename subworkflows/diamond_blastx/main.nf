process createDB {

    label 'diamond'

    tag "building ${id} database"

    input:
    tuple val(id), path(proteins_path)

    output:
    tuple val(id), path(diamond_db)

    script:
    diamond_db = "${id}.dmnd"
    """
        echo "Building ${diamond_db} database"
        diamond makedb --in ${proteins_path} -d ${id};
    """
}

//get hsp in gff3
process runDiamond {

    publishDir(params.OUTPUT, mode : 'copy')
    label "diamond"

    // show in the log which input file is analysed
    tag "running ${id} against ${dmnd_database_file}"

    input:
    tuple val(id), path(genome), path(dmnd_database_file)

    output:
    tuple val(id), path(output_name)

    script:
    output_name = "${id}.hsp.gff"
    tmp = "${id}.hsp.out"
    """
    fmt6_custom='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qstrand qframe'
    diamond blastx --db ${dmnd_database_file} \
                    --query ${genome} \
                    --max-target-seqs 0 \
                    --max-hsps 0  \
                    --outfmt \$fmt6_custom \
                    --evalue 0.0001 \
                    --block-size 0.8 \
                    --threads ${task.cpus} \
                    --out ${tmp}

    awk 'BEGIN{OFS="\t"}{if (\$13=="-"){frame=-(\$14);print \$1,"blastx","hsp",\$8,\$7,\$12,\$13,frame,\$2}else if(\$13=="+"){print \$1,"blastx","hsp",\$7,\$8,\$12,\$13,\$14,\$2}}' \
                ${tmp} | sort -k1,1 -k4,5n -k9,9 > ${output_name}
    """
}

workflow diamond_blastx {
    
    take:
    genomes
    uniref_proteins

    main:

    diamond_dbs = createDB(uniref_proteins)

    diamond_input = genomes.join(diamond_dbs)

    diamond_output = runDiamond(diamond_input)
     
    emit:
    diamond_output
}
