
# coding: utf-8

# In[1]:

with open("pth1.txt") as f: text = f.read()


# In[3]:

import pydna


# In[4]:

seqs = pydna.parse(text)


# In[5]:

len(seqs)


# In[6]:

from ypkpathway import PathWay


# In[7]:

x= PathWay(seqs)


# In[9]:

x.generate_files()


# In[11]:

x.tp_gene_tp


# In[14]:

a=x.tp_gene_tp[0]


# In[17]:

a.right_product()


# In[ ]:



