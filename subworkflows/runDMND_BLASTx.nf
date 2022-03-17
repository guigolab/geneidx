/*
*  Geneid module.
*/

// Parameter definitions
// params.CONTAINER = "quay.io/biocontainers/diamond:2.0.14--hdcc8f71_0"
params.CONTAINER = "quay.io/biocontainers/diamond:0.9.30--h56fc30b_0"
// params.OUTPUT = "protein_DBs"
// params.LABEL = ""


/*
 * Uncompressing if needed
 */
process UncompressFASTA {

    // where to store the results and in which way
    // publishDir(params.OUTPUT, pattern : '*.fa')

    // // indicates to use as a container the value indicated in the parameter
    // container params.CONTAINER

    // show in the log which input file is analysed
    tag "${ref_to_index}"

    input:
    file (ref_to_index)

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
    // perl -i -lane 'if (/^>/) { (\$id, \$chr)=\$_=~/^>([\\w|.]+)[\\s\\w]+, [\\w]+: (\\w+)/; print ">".\$chr} else {print}' ${main_genome_file}
    // perl -i -lane 'if (/^>/) { ($id, $chr)=$_=~/^>([\w|.]+)[\s\w]+, chromosome: (\w+)/; print ">".$chr} else {print}' ${main_genome_file}
}



process runDIAMOND_getHSPs {

    // where to store the results and in which way
    // publishDir(params.OUTPUT, pattern : '*.out')

    // indicates to use as container the value indicated in the parameter
    container params.CONTAINER

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${main_genome_name} against ${main_database_name}"

    input:
    path(dmnd_database_file)
    path(reference_genome_file)

    output:
    path ("${main_genome_name}.${main_database_name}.hsp.out")

    script:
    main_database_name = dmnd_database_file.BaseName
    main_genome_name = reference_genome_file.BaseName
    """
    fmt6_custom='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qstrand qframe'
    diamond blastx --db ${dmnd_database_file} \
                   --query ${reference_genome_file} \
                   --max-target-seqs 0 \
                   --max-hsps 0  \
                   --outfmt \$fmt6_custom \
                   --evalue 0.0001 \
                   --block-size 2.0 \
                   --threads 4 \
                   --out ${main_genome_name}.${main_database_name}.hsp.out
    rm ${reference_genome_file}
    """
                   // --max-hsps 0  \
    // mkdir -p tmp;
    // --tmpdir ./tmp/ \
    // https://www.metagenomics.wiki/tools/blast/evalue
}


process matches_to_GFF {

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy', pattern : '*.gff')

    // default container already has awk
    // container "quay.io/guigolab/geneid:1.4.5"

    tag "building ${main_matches_name} database"


    input:
    path(matches_file)

    output:
    path ("${main_matches_name}.gff")

    script:
    main_matches_name = matches_file.BaseName
    """
    awk 'BEGIN{OFS="\t"}{if (\$13=="-"){frame=-(\$14);print \$1,"blastx","hsp",\$8,\$7,\$12,\$13,frame,\$2}else if(\$13=="+"){print \$1,"blastx","hsp",\$7,\$8,\$12,\$13,\$14,\$2}}' ${main_matches_name}.out | sort -k1,1 -k4,5n -k9,9 | egrep -v '^MT' | cut -d '|' -f1 > ${main_matches_name}.gff
    """
}



/*
 * Workflow connecting the different pieces
 */

workflow alignGenome_Proteins {

    // definition of input
    take:
    prot_DB_file
    genome_file

    main:

    /*
    * Compressed file as input
    */
    genome_filename = UncompressFASTA(genome_file)

    matches = runDIAMOND_getHSPs(prot_DB_file, genome_filename)


    /*
    * Uncompressed file as input

    matches = runDIAMOND_getHSPs(prot_DB_file, genome_file)
    */

    // Merge the BLASTx matches with each other
    matchesGFF = matches_to_GFF(matches)

    emit:
    matchesGFF

}
