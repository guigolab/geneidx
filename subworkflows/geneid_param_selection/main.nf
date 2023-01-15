/*
*  Get parameters module.
*/

/*
 * Defining the output folders.
 */
OutputFolder = "${params.output}"

/*
 * Workflow for choosing and providing the parameter file
 */

process getMissingParams {
    input:
        val params
        val param_keys
    output:
        val missing_params
    exec:
        keys_from_user = params.collect{ it.key }
        missing_params = param_keys.findAll { !keys_from_user.contains(it) }
}

process getPathFromClosestTaxon {
    label 'geneidx'

    input:
        val target_lineage
        path matrix

    output:
        stdout emit: param_file
    
    script:
    params_path = params.auto_params_selection_files_path
    """
    #!/usr/bin/env python3
    # coding: utf-8
    import pandas as pd

    # Define an alternative in case everything fails
    selected_param = "Homo_sapiens.9606.param"


    data = pd.read_csv("${matrix}", sep = "\t")

    def split_n_convert(x):
        return [int(i) for i in x.replace("'", "").strip("[]").split(", ")]

    data.loc[:,"taxidlist"] = data.loc[:,"taxidlist"].apply(split_n_convert)

    ###
    # Separate the lineages into a single taxid per row
    ###
    exploded_df = data.explode("taxidlist")
    exploded_df.columns = ["species", "taxid_sp", "parameter_file", "taxid"]
    exploded_df.loc[:,"taxid"] = exploded_df.loc[:,"taxid"].astype(int)

    ###
    # Get the species of interest lineage
    ###
    query = pd.DataFrame(${target_lineage})
    query.columns = ["taxid"]
    query.loc[:,"taxid"] = query.loc[:,"taxid"].astype(int)

    ###
    # Intersect the species lineage with the dataframe of taxids for parameters
    ###
    intersected_params = query.merge(exploded_df, on = "taxid")

    ###
    # If there is an intersection, select the parameter whose taxid appears
    #   less times, less frequency implies more closeness
    ###
    if intersected_params.shape[0] > 0:

        taxid_closest_param = intersected_params.loc[:,"taxid"].value_counts().sort_values().index[0]

        selected_param = intersected_params[intersected_params["taxid"] == taxid_closest_param].loc[:,"parameter_file"].iloc[0]
    
    print("${params_path}", selected_param, sep = "/", end = "")
    """
}

process getLineage {
    input:
        val taxid
    output:
        val lineage
    exec:
        response = new URL("https://www.ebi.ac.uk/ena/browser/api/xml/${taxid}?download=false").text
        xml = new XmlSlurper().parseText(response)
        lineage = xml[0].children()[0].children()[0].children().collect { it.attributes()['taxId']}
}


workflow geneid_param_selection {

    // definition of input
    take:
    taxid

    main:

    lineage = getLineage(taxid) 
    matrix = file(params.auto_params_selection_matrix_path)
    selected_parameter_file = getPathFromClosestTaxon(lineage, matrix) 

    emit:
    param_file = selected_parameter_file

}
