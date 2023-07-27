

process getProtRepeats {

    // // Here we could add the option for a possible output
    // //   of the .tsv file
    // ${data_path}/

    // indicates to use as a container the value indicated in the parameter
    // container "geneidx"

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"

    // show in the log which input file is analysed
    tag "${taxid}"

    input:
    tuple val(id), val(taxid)
    val update // unless you want to look for new repeat files, leave it as False or 0
    path data_path

    output:
    tuple val(id), stdout


    script:
    """
    #!/usr/bin/env python3
    # coding: utf-8

    import os, sys
    import pandas as pd
    import requests
    from lxml import etree


    # Define an alternative in case everything fails
    #### MISSING ####


    # Define functions
    def choose_node(main_node, sp_name):
        for i in range(len(main_node)):
            if main_node[i].attrib["scientificName"] == sp_name:
                return main_node[i]
        return None


    # given a node labelled with the species, and with
    # lineage inside it returns the full path of the lineage
    def sp_to_lineage_clean ( sp_sel ):
        lineage = []

        if sp_sel is not None:
            lineage.append(sp_sel.attrib["taxId"])

        for taxon in sp_sel:
            if taxon.tag == 'lineage':
                lin_pos = 0
                for node in taxon:
                    if "rank" in node.attrib.keys():
                        lineage.append( node.attrib["taxId"] )
                    else:
                        lineage.append( node.attrib["taxId"] )
                    lin_pos += 1
        return(lineage)


    #
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


    # given a species from Ensembl, receive the link to the repeats file
    # as well as the filename alone
    def species_to_repeats(spec):

        species_specific_link = f"http://ftp.ebi.ac.uk/pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/{spec}/"
        species_specific_data = pd.read_html(species_specific_link)[0]

        species_specific_data = species_specific_data[species_specific_data["Name"].isna().astype(int) == 0]
        species_specific_data = species_specific_data[species_specific_data["Name"].str.endswith(".fa")
                                                     ].reset_index(drop = True)

        species_specific_data = species_specific_data[["Name", "Last modified", "Size"]]

        species_specific_data["Last modified"] = pd.to_datetime(species_specific_data["Last modified"])

        selected_repeats = species_specific_data.sort_values(by = ["Last modified", "Name"],
                                                         ascending = False).iloc[0]["Name"]

        return (species_specific_link, selected_repeats)



    if ${update}:
        ###
        # We want to update the lists as new repeat files may have appeared
        ###

        data_species = pd.read_html('http://ftp.ebi.ac.uk/pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/')[0]
        data_species = data_species[data_species["Name"].isna().astype(int) == 0]
        data_species = data_species[data_species["Name"].str.endswith("/")].reset_index(drop = True)

        # Select the list of species
        list_species = list(data_species["Name"].str.strip(" /").str.capitalize().values)

        # Put the list into a dataframe
        data_names = pd.DataFrame(list_species, columns = ["Species_name"])

        # Generate the dataframe with the filename and lineage information
        list_repeats_lineages = []
        for species_no_space in data_names.loc[:,"Species_name"] :
            species = species_no_space.replace("_", " ")

            # Get lineage
            lineage_sp = getLineage_from_Name(species)

            if lineage_sp is not None:
                list_repeats_lineages.append((species_no_space.lower(), lineage_sp))


        # Put the information into a dataframe
        data = pd.DataFrame(list_repeats_lineages, columns = ["species", "taxidlist"])

        data.to_csv("repeats_df.tsv", sep = "\t", index = False)
        # print("New parameters saved")


    else:
        ###
        # We want to load the previously generated dataframe
        ###
        data = pd.read_csv("${data_path}/repeats_df.tsv", sep = "\t")

        def split_n_convert(x):
            return [int(i) for i in x.replace("'", "").strip("[]").split(", ")]
        data.loc[:,"taxidlist"] = data.loc[:,"taxidlist"].apply(split_n_convert)

    # Following either one or the other strategy we now have N parameters to choose.
    # print(data.shape[0], "parameters available to choose")



    ###
    # Separate the lineages into a single taxid per row
    ###
    exploded_df = data.explode("taxidlist")
    exploded_df.columns = ["species", "taxid"]
    #exploded_df.loc[:,"taxid"] = exploded_df.loc[:,"taxid"].astype(int)

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

        selected_sp_repeats = intersected_params[intersected_params["taxid"] == taxid_closest_param].loc[:,"species"].iloc[0]

        link, file_sp = species_to_repeats(selected_sp_repeats)
        print(file_sp + "\\t" + link+file_sp, end = '')

    """

}





process downloadRepFasta {

    // where to store the results and in which way
    publishDir(params.OUTPUT, mode : 'copy', pattern : '*.fa')

    // indicates to use as a container the value indicated in the parameter
    // container "ferriolcalvet/geneidx"

    // indicates to use as a label the value indicated in the parameter
    label "geneidx"
    // show in the log which input file is analysed
    tag "${repeats_filename}"

    input:
    tuple val(id), val(repeats_desc)

    output:
    tuple val(id), path(repeats_filename)

    script:
    (repeats_file_desc, repeats_filename ) = (repeats_desc =~ /([A-Za-z\d\._\+]+)\t+(.*)/)[0]

    """
    #!/usr/bin/env python

    import requests, os

    info = "${repeats_desc}".strip().split("\\t")

    rep_filename = info[0]
    rep_filelink = info[1]


    if os.path.exists(f"${params.repeats_data_path}/{rep_filename}") :
      print(f"{rep_filename} already downloaded")
      os.symlink( f"${params.repeats_data_path}/{rep_filename}", rep_filename );

    else:
      print(f"Downloading {rep_filename}")
      r = requests.get(rep_filelink, allow_redirects=True)
      open(rep_filename, 'wb').write(r.content)
    """

}




/*
 * Workflow for obtaining the estimates of the exon sequences
 */

workflow repeat_download {

    // definition of input
    take:
    metadata
    data_path

    main:
    repeat_file_down = getProtRepeats(metadata, params.update_repeat_prot_file, data_path) | downloadRepFasta

    emit:
    repeat_file_down

}