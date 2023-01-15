import sys
import pandas as pd
import getopt

def parse_arguments(argv):
    arguments = dict()
    arg_help = "{0} -m <matrix> -l <lineage>".format(argv[0])

    try:
        opts, args = getopt.getopt(argv[1:], "l:m:p:", ["lineage=", 
        "matrix=","params_path="])
    except:
        print(arg_help)
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ("-l", "--lineage"):
            lineage_list = map(int, arg.strip('[]').split(','))
            arguments['lineage'] = lineage_list
        elif opt in ("-m", "--matrix"):
            arguments['matrix'] = arg
        elif opt in ("-p", "--params_path"):
            arguments['params_path'] = arg
    return arguments
# Define an alternative in case everything fails
selected_param = "Homo_sapiens.9606.param"

arguments = parse_arguments(sys.argv)

data = pd.read_csv(arguments['matrix'], sep = "\t")

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
query = pd.DataFrame(arguments['lineage'])
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

print(f"{arguments['params_path']}/{selected_param}")
