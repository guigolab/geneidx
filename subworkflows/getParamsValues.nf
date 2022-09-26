/*
*  Get parameters module.
*/

// Parameter definitions
params.CONTAINER = "ferriolcalvet/geneidx"
// params.OUTPUT = "geneid_output"
// params.LABEL = ""


/*
 * Defining the output folders.
 */
OutputFolder = "${params.output}"



process getParamName {

    // indicates to use as a container the value indicated in the parameter
    container "ferriolcalvet/geneidx"
    // container "ferriolcalvet/python-geneid-params"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${taxid}"

    input:
    val taxid
    val update // unless you manually add a new parameter file, leave it as False
    path param_path

    output:
    stdout emit: description


    script:
    """
    #!/usr/bin/env python3
    # coding: utf-8

    import os, sys
    import pandas as pd
    import requests
    from lxml import etree


    # Define an alternative in case everything fails
    selected_param = "Homo_sapiens.9606.param"


    # Define functions
    def choose_node(main_node, sp_name):
        for i in range(len(main_node)):
            if main_node[i].attrib["scientificName"] == sp_name:
                #print(main_node[i].attrib["rank"],
                #    main_node[i].attrib["taxId"],
                #    main_node[i].attrib["scientificName"])
                return main_node[i]
        return None


    # given a node labelled with the species, and with
    # lineage inside it returns the full path of the lineage
    def sp_to_lineage_clean ( sp_sel ):
        lineage = []

        if sp_sel is not None:
            lineage.append(sp_sel.attrib["taxId"])

        for taxon in sp_sel:
            #print(taxon.tag, taxon.attrib)
            if taxon.tag == 'lineage':
                lin_pos = 0
                for node in taxon:
                    if "rank" in node.attrib.keys():
                        lineage.append( node.attrib["taxId"] )
                    else:
                        lineage.append( node.attrib["taxId"] )
                    lin_pos += 1
        return(lineage)



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


    if ${update}:
        ###
        # We want to update the lists as new parameters may have been added
        ###

        # List files in directory
        list_species_taxid_params = os.listdir('${param_path}/*.param')
        list_species_taxid = [x.split('.')[:2] for x in list_species_taxid_params]


        # Put the list into a dataframe
        data_names = pd.DataFrame(list_species_taxid, columns = ["Species_name", "taxid"])


        # Generate the dataframe with the filename and lineage information
        list_repeats_taxids = []
        for species_no_space, taxid in zip(data_names.loc[:,"Species_name"], data_names.loc[:,"taxid"]):
            species = species_no_space.replace("_", " ")
            response = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/textsearch?domain=taxon&query={species}")
            xml = response.content
            if xml is None or len(xml) == 0:
                continue

            root = etree.fromstring(xml)
        #     print(species)
            sp_sel = choose_node(root, species)
            if sp_sel is None:
                continue
        #     print(sp_sel.attrib.items())getParamName
            lineage_sp = sp_to_lineage_clean(sp_sel)

            param_species = f"{species_no_space}.{taxid}.param"
            list_repeats_taxids.append((species, taxid, param_species, lineage_sp))
            # print((ens_sp, species, link_species, lineage_sp))


        # Put the information into a dataframe
        data = pd.DataFrame(list_repeats_taxids, columns = ["species", "taxid", "parameter_file", "taxidlist"])

        data.to_csv("/data/Parameter_files.taxid/params_df.tsv", sep = "\t", index = False)
        # print("New parameters saved")



    else:
        ###
        # We want to load the previously generated dataframe
        ###
        data = pd.read_csv("${param_path}/params_df.tsv", sep = "\t")

        def split_n_convert(x):
            return [int(i) for i in x.replace("'", "").strip("[]").split(", ")]
        data.loc[:,"taxidlist"] = data.loc[:,"taxidlist"].apply(split_n_convert)

    # Following either one or the other strategy we now have N parameters to choose.
    # print(data.shape[0], "parameters available to choose")



    ###
    # Separate the lineages into a single taxid per row
    ###
    exploded_df = data.explode("taxidlist")
    exploded_df.columns = ["species", "taxid_sp", "parameter_file", "taxid"]
    exploded_df.loc[:,"taxid"] = exploded_df.loc[:,"taxid"].astype(int)

    ###
    # Get the species of interest lineage
    ###
    query = pd.DataFrame(get_organism(int(${taxid})))
    query.columns = ["taxid"]
    query.loc[:,"taxid"] = query.loc[:,"taxid"].astype(int)
    # print(query)

    ###
    # Intersect the species lineage with the dataframe of taxids for parameters
    ###
    intersected_params = query.merge(exploded_df, on = "taxid")
    # print(intersected_params.shape)

    ###
    # If there is an intersection, select the parameter whose taxid appears
    #   less times, less frequency implies more closeness
    ###
    if intersected_params.shape[0] > 0:
        #print(intersected_params.loc[:,"taxid"].value_counts().sort_values())

        taxid_closest_param = intersected_params.loc[:,"taxid"].value_counts().sort_values().index[0]
        #print(taxid_closest_param)

        selected_param = intersected_params[intersected_params["taxid"] == taxid_closest_param].loc[:,"parameter_file"].iloc[0]
        print("${params.parameter_path}", selected_param, sep = "/", end = '')
    """

}






process paramSplitValues {

    // where to store the results and in which way
    // publishDir(params.OUTPUT, mode : 'copy', pattern : '*.gff3')

    // indicates to use as a container the value indicated in the parameter
    container "ferriolcalvet/geneidx"
    // container "ferriolcalvet/python-geneid-params"

    // indicates to use as a label the value indicated in the parameter
    label (params.LABEL)

    // show in the log which input file is analysed
    tag "${param_file_name}"

    input:
    path param_file
    val params_available

    output:
    stdout emit: list_params


    script:
    param_file_name = param_file.BaseName
    """
    #!/usr/bin/env python3
    # coding: utf-8

    list_all = set(['absolute_cutoff_exons', 'coding_cutoff_oligos', 'no_score',
                    'site_factor', 'exon_factor', 'hsp_factor', 'exon_weights'])

    dict_all_params_with_spaces = eval(\"\"\"${params_available}\"\"\".replace("[","{\\"").replace("]","}").replace(",",",\\"").replace(":","\\":"))
    dict_all_params = { key.replace(' ', ''): value for key, value in dict_all_params_with_spaces.items()}
    #print(dict_all_params)

    def find_missing_param(list_defined, param_file):

        params_to_search = list(list_all - set(list(list_defined)))

        for text_to_find in params_to_search:

            started = 0

            with open(param_file) as f:
                for line in f:
                    if started == 0:
                        if line.strip().lower() == text_to_find:
                            param_name = line.strip().lower()
                            #print(line.strip())
                            started = -1

                    # we have just found the label of the parameter of interest
                    # here we need read the value
                    elif started == -1:
                        vals = line.strip().split(' ')

                        # we check whether there are multiple values,
                        # if so, we take the value of the last one
                        if vals.count(vals[0]) > vals.count(vals[-1]):
                            dict_all_params[param_name] = float(vals[0])
                            #print(vals[0])

                        else:
                            dict_all_params[param_name] = float(vals[-1])
                            #print(vals[-1])

                        started = 0

                        break
        # return str(dict_all_params).replace("{", "[").replace("}", "]").replace("'", "\\"")
        return str(dict_all_params).replace("{", "").replace("}", "").replace("'", "").replace(" ", "")

    print(find_missing_param(list(dict_all_params.keys()), "${param_file}"), end = "")
    """
}


/*
 * Workflow for choosing and providing the parameter file
 */

workflow param_value_selection_workflow {
    // definition of input
    take:
    taxid
    update_list
    params_path
    param_list

    main:
    param_file_down = getParamName(taxid, update_list, params_path)

    // channel from list of params to find and report then here
    param_vals_out = paramSplitValues(param_file_down, param_list)



    emit:
    params_map = param_vals_out.list_params

}
