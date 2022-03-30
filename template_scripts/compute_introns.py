#!/usr/bin/env python
# coding: utf-8

# In[75]:


import sys
import pandas as pd


# """
# awk 'OFS="\t"{print $1, $4-40, $5+40, $9$7}' Homo_sapiens.GRCh38.dna_rm.primary_assembly.Vertebrates.5references.22.hsp.gff | sort -k1,1 -k4,4 -k2,2n > Homo_sapiens.GRCh38.dna_rm.primary_assembly.Vertebrates.5references.22.hsp.resorted.gff
# """

# In[76]:


#input_file = "Homo_sapiens.GRCh38.dna_rm.primary_assembly.Vertebrates.5references.22.hsp.resorted.gff"
#output_file = "Homo_sapiens.GRCh38.dna_rm.primary_assembly.Vertebrates.5references.22.hsp.introns.gff3"
#max_intron_size = 10000
input_file = sys.argv[1]
output_file = sys.argv[2]
max_intron_size = int(sys.argv[3])


# In[77]:


data = pd.read_csv(input_file,
                   sep = "\t", header = None)


# In[78]:


#print(data.shape)
#data.head(10)


# In[79]:


prev_id = ""
prev_end = 0
intron_l = []

for index, row in data.iterrows():
#     dif = row[2] - prev_end
    if row[3] == prev_id:
        intron_l.append(list(row) + [prev_end])
    else:
        prev_id = row[3]
    prev_end = row[2]


# In[80]:


dd_intron = pd.DataFrame(intron_l)
dd_intron = dd_intron[[0,4,1,3]]
dd_intron[3] = dd_intron[3].apply(lambda x : x[-1])
dd_intron = dd_intron[(dd_intron[1] - dd_intron[4]) > 0]
dd_intron = dd_intron[(dd_intron[1] - dd_intron[4]) < max_intron_size].reset_index(drop = True).reset_index()
#print(dd_intron.shape)
#dd_intron.head()


# In[81]:


dd_intron["dot"] = '.'
dd_intron["source"] = 'hsp'
dd_intron["type"] = 'CDS'


# In[82]:


dd_intron["index"] = "ID=" + dd_intron["index"].astype(str) + ";Parent=" + dd_intron["index"].astype(str) + ";"
dd_intron = dd_intron[[0, 'source', 'type', 4, 1, 'dot', 3, 'dot', 'index']]
#dd_intron.head()


# In[84]:


dd_intron.to_csv(output_file,
                 sep = "\t",
                 header = None,
                 index = None)


# awk '!found[$1"\t"$2"\t"$3"\t"$5]++'
