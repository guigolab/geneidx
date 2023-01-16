// Parameter definitions
// params.CONTAINER = "quay.io/biocontainers/diamond:2.0.14--hdcc8f71_0"


// /*
//  * Defining the output folders.
//  */
// OutputFolder = "${params.output}"
// OutputFolderInternal = "${OutputFolder}/internal"


/*
 * Defining the module / subworkflow path, and include the elements
 */
subwork_folder = "${projectDir}/subworkflows/"

include { unzipFasta } from "${subwork_folder}/tools"


process createDB {

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy', pattern : '*.dmnd')

    // indicates to use as a label the value indicated in the parameter
    label 'diamond'

    // show in the log which input file is analysed
    // tag "${ref}"
    tag "building ${main_proteins_name} database"

    input:
    path(reference_proteins_file)

    output:
    path ("${main_proteins_name}.dmnd")

    script:
    main_proteins_name = reference_proteins_file.BaseName
    // we used this before when we were not cleaning the fasta identifiers
    // query_curated = query.toString().tokenize('|').get(1)
    """
    echo "${params.OUTPUT}/${main_proteins_name}.dmnd"
    if [ ! -s ${params.OUTPUT}/${main_proteins_name}.dmnd ]; then
        echo "Building ${main_proteins_name}.dmnd database"
        diamond makedb --in ${reference_proteins_file} -d ${main_proteins_name};
    else
        echo "${main_proteins_name}.dmnd database already built"
        ln -s ${params.OUTPUT}/${main_proteins_name}.dmnd ${main_proteins_name}.dmnd;
    fi
    rm ${reference_proteins_file};
    """
}



/*
 * Workflow connecting the different pieces
 */

workflow diamond_db_build {
    take:
    prot_file

    main:
    proteins_filename = unzipFasta(prot_file)

    protein_DB = createDB(proteins_filename)

    emit:
    protein_DB
}
