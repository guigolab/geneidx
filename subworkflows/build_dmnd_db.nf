/*
*  Geneid module.
*/

// Parameter definitions
// params.CONTAINER = "geneid_path"
params.CONTAINER = "quay.io/biocontainers/bowtie:1.2.3--py37hc9558a2_0"
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



// /*
//  * Indexing if needed
//  */
//
// process Index {
//
//     // where to store the results and in which way
//     // publishDir(params.OUTPUT, pattern : '*.fa.i')
//
//     // // indicates to use as a container the value indicated in the parameter
//     // container params.CONTAINER
//
//     // show in the log which input file is analysed
//     tag "${main_genome_file}"
//
//     input:
//     path main_genome_file
//
//     output:
//     path ("${main_genome_file}.i")
//
//     script:
//     """
//     if [ ! -s  ${main_genome_file}.i ]; then
//         echo "indexing genome ${main_genome_file}"
//         fastaindex -f ${main_genome_file} -i ${main_genome_file}.i
//     fi
//     """
//     // cut -d ' ' -f1 ${main_genome_file}.i >> ${main_genome_file}.list
// }



process runDIAMOND_makedb {

    // where to store the results and in which way
    publishDir(params.OUTPUT, pattern : '*.gff3')

    // indicates to use as container the value indicated in the parameter
    // container params.CONTAINER

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
    diamond makedb --in ${reference_proteins_file} -d $prot_dbs/${main_proteins_name}
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

    // main part where we connect two modules, indexing and predicting
    main:


    // I SHOULD ADD A CONDITION HERE TO CHECK IF THE DMND FILE
    //    EXISTS IN THE CORRESPONDING FOLDER

    proteins_filename = UncompressFASTA(prot_file)
    // genome_filename.subscribe {  println "Got: $it"  }

    // index_filename = Index(genome_filename)
    // // index_filename.subscribe {  println "Got: $it"  }

    // geneid_param.view()
    // genome_filename.view()
    // index_filename.view()

    // we call the runGeneid_fetching module using the channel for the queries
    protein_DB = runDIAMOND_makedb(proteins_filename)


    emit:
    // index_filename
    protein_DB
}
