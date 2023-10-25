import groovy.json.JsonSlurper

process getParamValues {

    // indicates to use as a label the value indicated in the parameter
    label 'geneidx'

    // show in the log which input file is analysed
    tag "${param_file_name}"

    input:
    tuple val(id), val(taxid), path(param_file)

    output:
    tuple val(id), stdout

    script:
    params_available = "${params.maps_param_values}"
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

process getFileFromTaxon {
    label 'geneidx'

    input:
        tuple val(id), val(taxid), val(target_lineage)
        path matrix

    output:
        tuple val(id), val(taxid), stdout
    
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
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(id), val(taxid)
    output:
        tuple val(id), val(taxid), val(lineage)
    exec:
        response = new URL("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/${taxid}").text
        json = new JsonSlurper().parseText(response)
        lineage = json.taxonomy_nodes[0].taxonomy.lineage.reverse()
        return lineage
}


process splitParams {

    // indicates to use as a label the value indicated in the parameter
    label 'geneidx'

    input:
    tuple val(id), val(taxid), path(param_file)

    output:
    tuple val(id), val(taxid), path("${param_file_name}.acceptor_profile.param"), path ("${param_file_name}.donor_profile.param"), path ("${param_file_name}.start_profile.param"), path ("${param_file_name}.stop_profile.param")

    script:
    param_file_name = param_file.getName()
    """
    #!/usr/bin/env python3
    # coding: utf-8

    import re

    pattern = r'([a-zA-Z]+)_[Pp]rofile'

    started = 0
    files_created = []

    # removing the new line characters
    with open("${param_file_name}") as f:
        for line in f:
            if started == 0:
                matches = re.match(pattern, line)
                # if there is any match, get the groups and start a new file
                if matches is not None:
                    groups = matches.groups()
                    # create the file that will contain the given profile
                    # only if it has not been created before
                    filename = "${param_file_name}." + groups[0].lower() + "_profile.param"
                    if filename not in files_created:
                        files_created.append(filename)
                        started = -1
                        fW = open(filename, "w")
                        fW.write(line)

            # we have just found the label of the profile of interest
            # here we need to know the length and order
            elif started == -1:
                # write line to file
                fW.write(line)

                # process the content of this line
                clean_line = line.strip()
                splitted_line = clean_line.split(" ")
                length_pwm = int(splitted_line[0])
                order_pwm = int(splitted_line[3])

                rep_position = 4 ** (order_pwm + 1)
                total_to_read = rep_position * length_pwm

                # define started as the total number of lines to be read and written to the output file
                started = total_to_read + 1

            # these are the lines that we are reading and adding to the output file of the profile
            elif started > 1:
                # write line to file and substract one
                fW.write(line)
                started -= 1
                # print(line, end = '')

            # this is the last line to be added, so we add it and close the file
            elif started == 1:
                # write line to file and substract one
                fW.write(line)
                started -= 1
                # print(line, end = '')

                # close the file as we finished adding that profile
                fW.close()
    """
}

workflow geneid_param_creation {

    take:
    meta

    main:
    
    lineage = getLineage(meta)

    matrix = file(params.auto_params_selection_matrix_path)

    selected_parameter_file = getFileFromTaxon(lineage, matrix)

    splitted_params = splitParams(selected_parameter_file)

    param_values = getParamValues(selected_parameter_file)

    //emit tuple val(id), path(acceptor), path(donor), path(start), path(stop), stdout params
    geneid_parameters = splitted_params.join(param_values)
    
    emit:
    geneid_parameters
    // acc_pwm
    // don_pwm
    // sta_pwm
    // sto_pwm
    // param_valuesgetParamValues
    // id
}
