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

process DownloadFASTA_fromID {

    // where to store the results and in which way
    // publishDir(params.OUTPUT, pattern : '*.fa.gz')

    // // indicates to use as a container the value indicated in the parameter
    // container params.CONTAINER

    // show in the log which input file is analysed
    tag "${assembly} from ${name}"

    input:
    tuple val(assembly),val(chromosomes),val(taxid),val(name),path(param)

    output:
    path("${species_name}.${assembly}.fa.gz"), emit: genome
    path(param), emit: param_f

    script:
    species_name = name.replaceAll(" ", "_")
    link_to_ref = chromosomes.join(",")
    """
    wget "https://www.ebi.ac.uk/ena/browser/api/fasta/${link_to_ref}?download=true&gzip=true" -O ${species_name}.${assembly}.fa.gz;
    """
}
