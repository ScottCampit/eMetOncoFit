"""
Proteomics data analysis

- Read in mass spec data
- Determine methylation and acetylation quantities

@author: Scott Campit
"""

import os
import re
import pandas as pd

# Read in Krishna's protein expression data and figure out how he got these numbers later.
path = r"./../raw/PROTEOME/"
folder = os.listdir(path)
for fil in folder:
    tmp = pd.read_csv(path+fil, sep='\t')
    id, _ = os.path.splitext(fil)
    tissue = id.split('_')[0]
    cell = id.split('_')[-1]

    tmp.columns = ['Genes', 'Value']
    tmp = tmp.fillna(method='ffill')

    # restack genes
    tmp2 = pd.DataFrame(tmp['Genes'].str.split(';').tolist(), index=tmp['Value']).stack()
    tmp2 = tmp2.reset_index()[[0, 'Value']]
    tmp2.columns = ['Genes', 'Value']

    proteomics = tmp2.groupby('Genes')['Value'].median()
    proteomics = proteomics.reset_index()
    proteomics['Cancer'] = tissue
    proteomics['Cell Line'] = cell
    print(proteomics)


#import padua

#df = padua.io.read_maxquant(r'./prot/proteinGroups.txt')
#df = padua.filters.remove_reverse(df)
#df = padua.filters.remove_potential_contaminants(df)
#df = df.set_index('Protein', inplace=True)
#df = df.filter(regex='Intensity')
#df = padua.process.expand_side_table(df)
#df = df.filter(regex='Intensity ')
#print(df)
