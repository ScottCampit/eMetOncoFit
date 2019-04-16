"""
Proteomics data analysis

- Read in mass spec data
- Determine methylation and acetylation quantities

"""

import padua

df = padua.io.read_maxquant(r'./prot/proteinGroups.txt')
#df = padua.filters.remove_reverse(df)
#df = padua.filters.remove_potential_contaminants(df)
#df = df.set_index('Protein', inplace=True)
#df = df.filter(regex='Intensity')
#df = padua.process.expand_side_table(df)
#df = df.filter(regex='Intensity ')
print(df)
