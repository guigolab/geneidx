// Parameter definitions
// params.CONTAINER = "geneid_path"
params.OUTPUT = "geneid_output"
params.LABEL = ""

/*
 * choose the parameter file among the available ones
 */

process chooseParamFile {

    // container 'python:3.6-stretch'

    // where to store the results and in which way
    publishDir(params.OUTPUT)

    // // indicates to use as a container the value indicated in the parameter
    // container params.CONTAINER

    // show in the log which input file is analysed
    tag "Checking most suitable parameter file"

    input:
    val query_taxid

    output:
    // val("${link_to_ref}"), emit: fasta_name
    stdout emit: param_file

    script:
    """
    #!/usr/bin/env python
    import requests
    from lxml import etree

    import pandas as pd

    taxid_of_interest = $query_taxid
    table_with_param_data = "$projectDir/data/select-param-files/taxid_to_links.tsv"

    # we should remove this line at some point, but by now I have it to be safe
    selected_param_sp = "Homo_sapiens"
    selected_param = "https://genome.crg.es/pub/software/geneid/human.070123.param"
    selected_param = 9606


    # table_with_param_data = "$projectDir/data/select-param-files/taxid_to_links.tsv"
    scored_links = pd.read_csv(table_with_param_data, header = True, sep = "\t")


    # we use this function developed by Emilio R that takes advantage of the ENA API
    # to get the taxid lineage of an organism given a taxid
    def get_organism(taxon_id):
        response = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{taxon_id}?download=false") ##
        if response.status_code == 200:
            root = etree.fromstring(response.content)
            species = root[0].attrib
            lineage = []
            for taxon in root[0]:
                if taxon.tag == 'lineage':
                    for node in taxon:
                        lineage.append(node.attrib["taxId"])
        return lineage



    query = pd.DataFrame(get_organism(taxid_of_interest))
    query.columns = ["taxid"]
    query.loc[:,"taxid"] = query.loc[:,"taxid"].astype(int)
    intersected_params = query.merge(scored_links, on = "taxid").sort_values(by = "rank_pos", ascending = False)
    if intersected_params.shape[0] > 0:
        selected_param = intersected_params.loc[0,"param_taxid"]
        #selected_param = intersected_params.loc[0,"link"]
        selected_param_sp = intersected_params.loc[0,"species"]

    print(selected_param)
    """
}
