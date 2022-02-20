/*
*  Geneid module.
*/

// Parameter definitions
params.CONTAINER = "quay.io/guigolab/geneid:1.4.5"
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
    """
    // perl -i -lane 'if (/^>/) { (\$id, \$chr)=\$_=~/^>([\\w|.]+)[\\s\\w]+, [\\w]+: (\\w+)/; print ">".\$chr} else {print}' ${main_genome_file}
    // perl -i -lane 'if (/^>/) { ($id, $chr)=$_=~/^>([\w|.]+)[\s\w]+, chromosome: (\w+)/; print ">".$chr} else {print}' ${main_genome_file}
}



/*
 * Indexing if needed
 */

process Index {

    // where to store the results and in which way
    // publishDir(params.OUTPUT, pattern : '*.fa.i')

    // // indicates to use as a container the value indicated in the parameter
    // container params.CONTAINER

    // show in the log which input file is analysed
    tag "${main_genome_file}"

    input:
    path main_genome_file

    output:
    path ("${main_genome_file}.i")

    script:
    """
    if [ ! -s  ${main_genome_file}.i ]; then
        echo "indexing genome ${main_genome_file}"
        fastaindex -f ${main_genome_file} -i ${main_genome_file}.i
    fi
    """
    // cut -d ' ' -f1 ${main_genome_file}.i >> ${main_genome_file}.list
}



process runGeneid_fetching {

    // where to store the results and in which way
    publishDir(params.OUTPUT, pattern : '*.gff3')

    // indicates to use as container the value indicated in the parameter
    container params.CONTAINER

    // show in the log which input file is analysed
    // tag "${ref}"
    tag "run Geneid ${query}"

    // MAYBE WE CAN ADD SOMETHING ABOUT THE TASK.CPUS HERE ??

    input:
    path(reference_genome_file)
    path(reference_genome_index)
    path(geneid_param)
    val query

    output:
    path ("${main_genome_file}.*.gff3")

    script:
    main_genome_file = reference_genome_file.BaseName
    query_curated = query
    // we used this before when we were not cleaning the fasta identifiers
    // query_curated = query.toString().tokenize('|').get(1)
    """
    fastafetch -f ${reference_genome_file} -i ${reference_genome_index} -q \"${query}\" > ${main_genome_file}.${query_curated}
    geneid -3P ${geneid_param} ${main_genome_file}.${query_curated} > ${main_genome_file}.${query_curated}.gff3
    rm ${main_genome_file}.${query_curated}
    """
}



/*
 * Workflow connecting the different pieces
 */

workflow geneid_WORKFLOW {

    // definition of input
    take:
    ref_file
    geneid_param

    // main part where we connect two modules, indexing and predicting
    main:

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


    // ch -> 39 -> 40 -> mitochondrion
    // only process1 reads
    // ch -> mitochondrion
    // process 2 finishes wiht the download
    // ch -> mitochondrion -> 33 -> W
    // process 1 and 2 read from here


    ch.view()
    ch.subscribe {  println "Got: $it"  }

    // geneid_param.view()
    // genome_filename.view()
    // index_filename.view()

    // we call the runGeneid_fetching module using the channel for the queries
    // predictions = runGeneid_fetching(genome_filename.first(), index_filename.first(), geneid_param.first(), ch)


    // emit:
    // index_filename
    // pred = predictions
}
