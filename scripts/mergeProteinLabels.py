#!/usr/bin/env python
# coding: utf-8

# In[12]:


import pandas as pd
import sys


# In[13]:


data = pd.read_table(sys.argv[1],
                     names=["seq", "source", "feature", "start", "end", "strand",
                           "score", "phase", "attributes", "protein"],
                     header = None)
data.head()


# In[16]:


grouped_data = data.groupby( ["seq", "source", "feature", "start", "end",
                              "strand", "score", "phase",
                              "attributes"]).agg({"protein" : ','.join }).reset_index()


# In[17]:


#grouped_data.head()


# In[21]:


updated_attributes = [ x if y == '.' else x + f";ProteinIDs={y}" for x, y in zip(grouped_data["attributes"].values, grouped_data["protein"].values) ]


# In[24]:


grouped_data["attributes"] = updated_attributes
grouped_data.drop("protein", axis = 1, inplace = True)



grouped_data.sort_values(["seq", "start", "end"],
              ascending = [True, True, False], inplace = True)
# In[25]:


grouped_data.to_csv(sys.argv[2],
                       sep = "\t",
                       header = None,
                       index = None)


# In[ ]:
