process getParamValues {

    // indicates to use as a container the value indicated in the parameter
    // container "ferriolcalvet/python-geneid-params"

    // indicates to use as a label the value indicated in the parameter
    label 'geneidx'

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
workflow geneid_param_values {
    take:
    param_file
    param_list

    main:
    param_vals_out = getParamValues(param_file, param_list)

    emit:
    params_map = param_vals_out.list_params

}
