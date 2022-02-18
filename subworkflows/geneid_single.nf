/*
*  Geneid module.
*/

// Parameter definitions
// params.CONTAINER = "geneid_path"
// params.OUTPUT = "geneid_output"
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
    perl -i -lane 'if (/^>/) { (\$id, \$chr)=\$_=~/^>([\\w|.]+)[\\s\\w]+, [\\w]+: (\\w+)/; print ">".\$chr} else {print}' ${main_genome_file}
    """
    // perl -i -lane 'if (/^>/) { ($id, $chr)=$_=~/^>([\w|.]+)[\s\w]+, chromosome: (\w+)/; print ">".$chr} else {print}' ${main_genome_file}
}


process runGeneid_single {

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy', pattern : '*.gff3')

    // indicates to use as container the value indicated in the parameter
    // container params.CONTAINER

    // show in the log which input file is analysed
    // tag "${ref}"
    tag "run Geneid ${reference_genome_file}"

    // MAYBE WE CAN ADD SOMETHING ABOUT THE TASK.CPUS HERE ??

    input:
    path(reference_genome_file)
    path(geneid_param)

    output:
    path ("${main_genome_file}.gff3")

    script:
    main_genome_file = reference_genome_file.BaseName
    """
    geneid -3P ${geneid_param} ${reference_genome_file} > ${main_genome_file}.gff3
    rm ${reference_genome_file}
    """
}




/*
 * Workflow connecting the different pieces
 */

workflow geneid_WORKFLOW_single {

    // definition of input
    take:
    ref_file
    geneid_param

    // main part where we connect two modules, indexing and predicting
    main:

    genome_filename = UncompressFASTA(ref_file)
    // genome_filename.subscribe {  println "Got: $it"  }

    // we call the runGeneid_fetching module using the channel for the queries
    predictions = runGeneid_single(genome_filename, geneid_param)


    emit:
    pred = predictions
}
