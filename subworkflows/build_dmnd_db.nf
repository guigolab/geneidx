// Parameter definitions
params.CONTAINER = "quay.io/biocontainers/diamond:2.0.14--hdcc8f71_0"
params.OUTPUT = "protein_DBs"
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



process runDIAMOND_makedb {

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy', pattern : '*.dmnd')

    // indicates to use as container the value indicated in the parameter
    container params.CONTAINER

    // show in the log which input file is analysed
    // tag "${ref}"
    tag "building ${main_proteins_name} database"

    // MAYBE WE CAN ADD SOMETHING ABOUT THE TASK.CPUS HERE ??

    input:
    path(reference_proteins_file)

    output:
    path ("${main_proteins_name}.dmnd")

    script:
    main_proteins_name = reference_proteins_file.BaseName
    // we used this before when we were not cleaning the fasta identifiers
    // query_curated = query.toString().tokenize('|').get(1)
    """
    diamond makedb --in ${reference_proteins_file} -d ${main_proteins_name}
    rm ${reference_proteins_file}
    """
}



/*
 * Workflow connecting the different pieces
 */

workflow build_protein_DB {

    // definition of input
    take:
    prot_file

    main:
    // I SHOULD ADD A CONDITION HERE TO CHECK IF THE DMND FILE
    //    EXISTS IN THE CORRESPONDING FOLDER

    proteins_filename = UncompressFASTA(prot_file)
    // genome_filename.subscribe {  println "Got: $it"  }

    // we call the runGeneid_fetching module using the channel for the queries
    protein_DB = runDIAMOND_makedb(proteins_filename)


    emit:
    // index_filename
    protein_DB
}
