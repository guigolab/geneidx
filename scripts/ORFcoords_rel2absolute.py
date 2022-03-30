#!/usr/bin/python
import sys
import pandas as pd

orf_coords = sys.argv[1] # Relative ORF coordinates
gff3_file = sys.argv[2] # initial gff3
output = sys.argv[3] # output orf gff3

read_orf_coords = pd.read_csv(orf_coords, sep='\t', header=None)
read_gff3_file = pd.read_csv(gff3_file, sep='\t', header=None)

read_orf_coords.columns = ["id", "rel_start", "rel_end", "info", "frame2", "strand_ORF"]
read_gff3_file.columns = ["chr", "program", "region", "start", "end", "value", "strand", "frame", "id"]
read_orf_coords.loc[:,"id"] = 'ID='+ read_orf_coords.loc[:,"id"].astype(str) + ';Parent='+ read_orf_coords.loc[:,"id"].astype(str) + ';'

merged_orf_gff3 = read_gff3_file.merge(read_orf_coords, on='id', how='inner')




def fix_coords(x):
    if x["strand"] == '+':
        return ( int(x["start"] + x["rel_start"]), int(x["start"] + x["rel_end"]) - 1 )

    return ( int(x["end"] - x ["rel_end"])+ 1, int(x["end"] - x["rel_start"]) )

fixed_coords = merged_orf_gff3.apply(fix_coords, axis=1)
fixed_coords = pd.DataFrame(fixed_coords.to_list())
# print(fixed_coords)

merged_orf_gff3.loc[:,"start"] = fixed_coords.iloc[:,0]
merged_orf_gff3.loc[:,"end"]   = fixed_coords.iloc[:,1]

# print(merged_orf_gff3.columns)
ORF_coords = merged_orf_gff3[['chr', 'program', 'region', 'start', 'end', 'value', 'strand', 'frame', 'id']]
ORF_coords.to_csv(output, sep='\t', index=False, header=False)
# print(ORF_coords)
