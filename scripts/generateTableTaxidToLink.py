import numpy as np
import pandas as pd

taxa_ranks_list = ["Kingdom", "Subkingdom", "Phylum", "Subphylum", \
                    "Superclass", "Class", "Subclass", "Infraclass",  \
                    "Cohort", "Superorder", "Order", "Suborder", "Superfamily", \
                    "Family", "Subfamily", "Tribe", "Genus", "Subgenus", "Species", "Subspecies"]

data = pd.read_csv("../data/select-param-files/Link_Species_Taxids.txt", sep = "\t", header = None,
                   names = ["link", "species", "taxidlist"])
data.loc[:,"taxidlist"] = data.loc[:,"taxidlist"].str.strip().str.split(" ")


nodes_data = pd.read_csv("../data/select-param-files/nodes_interest.dmp", sep = "\t", header = None,
                         names = ["taxid", "parent_taxid", "rank"])
nodes_data.drop("parent_taxid", axis = 1, inplace = True)

# In[8]: DIFFERENT LISTS OF RANKS


ranks_list = ["superkingdom", "kingdom", "subkingdom", "phylum", "subphylum", "superclass", "class", "subclass",
              "infraclass", "cohort", "superorder", "order", "suborder", 'infraorder', 'parvorder', "superfamily",
              "family", "subfamily", "tribe", "genus", "subgenus", "species", "subspecies"]

# ranks_list2 = ['superkingdom', 'kingdom', 'phylum', 'subphylum', 'superclass', 'class',
#                'superorder', 'order', 'suborder', 'infraorder', 'parvorder', 'superfamily',
#                'family', 'subfamily', 'genus']
# ranks_to_keep = ["subphylum", "superclass", "class", "subclass", "infraclass", "cohort", "superorder",
#                  "order", "suborder", "superfamily", "family", "subfamily", "tribe", "genus",
#                  "subgenus", "species", "subspecies"]

# ranks_to_keep = ['superkingdom','kingdom','phylum','subphylum','class','order','family','genus','species']
# ranks_to_keep = ['subphylum','class','order','family','genus','species']
ranks_to_keep = ["subphylum", "superclass", "class", "subclass", "infraclass", "cohort",
                 "superorder", "order", "suborder", 'infraorder', 'parvorder', "superfamily",
                 "family", "subfamily", "tribe", "genus", "subgenus", "species", "subspecies"]

# we put the sorted ranks into a dataframe and we assign a value to each of them
# the subspecies rank has the higher value so the higher the value
# the more restricted the rank
rank_dat = pd.DataFrame(ranks_list).reset_index()
rank_dat.columns = ["rank_pos", "rank"]


selection_crit = [ x in ranks_to_keep for x in nodes_data.loc[:,"rank"] ]
nod_data = nodes_data[selection_crit].reset_index(drop = True)
# print(nod_data.shape)


def check_if_keep(dd):
    ddd = list(filter(lambda x: x != "", dd))
    sel_crit = [int(y) in nod_data.loc[:,"taxid"].values for y in ddd]
    return np.array(ddd)[sel_crit]

initial_taxid = [ x[0] for x in data.loc[:,"taxidlist"] ]
data.loc[:,"param_taxid"] = initial_taxid

filtered_taxid_list = data.loc[:,"taxidlist"].apply(check_if_keep)
data.loc[:,"taxidlist"] = filtered_taxid_list


# we explode the dataframe so that we have one taxid per line
exploded_df = data.explode("taxidlist")
exploded_df.columns = ["link", "species", "taxid", "param_taxid"]
exploded_df.loc[:,"taxid"] = exploded_df.loc[:,"taxid"].astype(int)

# we add the rank data to the data with links and lineage
annotated_links = exploded_df.merge(nod_data, on = "taxid")

# we add the score of the rank to the previous dataframe
scored_links = annotated_links.merge(rank_dat, on = "rank")


scored_links.to_csv("../data/select-param-files/taxid_to_links.tsv", header = True, sep = "\t", index = False)
