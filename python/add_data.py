# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 19:08:52 2018

Currently adding data to the current metabolic models:
- Km
- Metabolite concentrations
- miRNA
- lncRNA
- Thermodynamic information

@author: Scott Campit
"""

import sys, os, re, difflib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def metabolizer(input_path, output_path):
    '''

    '''

    # Import whole human km dataset
    os.chdir('C:\\Users\\scott\\Google Drive\\Work\\UM\\Research\\Sriram\\Data\\sabiodb\\sabioRK_Mammals')
    Human = pd.read_table('homosapien.txt')

    Human_uniprot = Human['UniprotID'].drop_duplicates(keep='first')
    np.savetxt(r'C:\Users\scott\humanUniprotID.txt', Human_uniprot, fmt='%s')
    genesymbol = pd.read_csv('C:\\Users\\scott\\Google Drive\\Work\\UM\\Research\\Sriram\\Data\\sabiodb\\sabioRK_Mammals\Human_GeneSymbol_Map.tab', sep='\t')
    Human = pd.merge(Human, genesymbol, how='left', left_on=Human.UniprotID, right_on=genesymbol.Entry)
    Human = Human.drop(columns=['Entry', 'Status'], axis=1)
    Human = Human[Human['EnzymeType'].str.contains('wildtype')]

    # Sabio RK dataset now only contains the WT
    Sabio_RK = Human[Human['parameter.type'].isin(['Km'])]
    Sabio_RK['log(Km)'] = np.log(Human['parameter.startValue'])
    Sabio_RK = Sabio_RK[np.isfinite(Sabio_RK['log(Km)'])]
    Sabio_RK['parameter.associatedSpecies'] = Sabio_RK['parameter.associatedSpecies'].str.lower()
    Sabio_RK['Substrate'] = Sabio_RK['Substrate'].str.lower()
    Sabio_RK['Product'] = Sabio_RK['Product'].str.lower()
    del Human, genesymbol, Human_uniprot

    # Separate GENES ----------------------------------------------------------

    # separate list of genes by spaces
    s = Sabio_RK['Gene names'].str.split(' ').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'Gene names'
    del Sabio_RK['Gene names']
    # drop duplicates
    Sabio_RK = Sabio_RK.join(s)
    Sabio_RK = Sabio_RK.drop_duplicates(keep='first')

    # separate list of genes by the forward slash
    s = Sabio_RK['Gene names'].str.split('/').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'Gene names'
    del Sabio_RK['Gene names']
    # drop duplicates
    Sabio_RK = Sabio_RK.join(s)
    Sabio_RK = Sabio_RK.drop_duplicates(keep='first')

    # Separate Substrates -----------------------------------------------------
    # separate list of genes by spaces
    s = Sabio_RK['Substrate'].str.split(';').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'Substrate'
    del Sabio_RK['Substrate']
    # drop duplicates
    Sabio_RK = Sabio_RK.join(s)
    Sabio_RK = Sabio_RK.drop_duplicates(keep='first')

    # Separate Products -------------------------------------------------------
    # separate list of genes by spaces
    s = Sabio_RK['Product'].str.split(';').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'Product'
    del Sabio_RK['Product']
    # drop duplicates
    Sabio_RK = Sabio_RK.join(s)
    Sabio_RK = Sabio_RK.drop_duplicates(keep='first')
    del s

    # get number of genes in Sabio RK dataset
    genes = Sabio_RK['Gene names'].unique()
    genes = pd.DataFrame(genes) # number is 1204

    # Read in NCI60 cell line
    os.chdir('C:\\Users\\scott\\Google Drive\\Work\\UM\\Research\\Sriram\\Data\\NCI60\\')
    nci_60 = pd.read_csv('WEB_DATA_METABOLON.txt')
    nci_60 = nci_60[nci_60['TITLE'].str.contains("X-") == False] # remove unannotated metabolites
    nci_60['TITLE'] = nci_60['TITLE'].str.lower()
    nci_60 = nci_60.set_index('TITLE')

    # separate cell lines with the '/'
    s = nci_60['cellname'].str.split('/').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'cellname'
    del nci_60['cellname']
    # drop duplicates
    nci_60 = nci_60.join(s)
    nci_60 = nci_60.drop_duplicates(keep='first')

    # extract important information
    clean_human_km = pd.DataFrame()
    clean_human_km['Gene'] = Sabio_RK['Gene names']
    clean_human_km['Metabolite'] = Sabio_RK['parameter.associatedSpecies']
    clean_human_km['log10(Km)'] = Sabio_RK['log(Km)']
    clean_human_km['Substrate'] = Sabio_RK['Substrate']
    clean_human_km['Product'] = Sabio_RK['Product']
    clean_human_km = clean_human_km.drop_duplicates(keep='first').dropna().fillna(0)
    clean_human_km = clean_human_km.dropna(axis=0)
    clean_human_km['log10(Km)'] = clean_human_km['log10(Km)'].fillna(0)

    # lower precision
    # multiplied by 100 and removed decimal
    clean_human_km['log10(Km)'] = clean_human_km['log10(Km)']*100
    clean_human_km = clean_human_km.astype({'log10(Km)': int})
    # multiplied by 100 and removed decimal
    nci_60 = nci_60[['VALUE', 'cellname']]
    nci_60.VALUE = nci_60.VALUE*100
    nci_60 = nci_60.astype({'VALUE': int})

    # Get metabolite values that are associated with biosynthesis -------------
    temp = pd.merge(nci_60, clean_human_km, how='inner', left_index=True, \
                    right_on='Product')

    # temp 1 contains log10(Km)
    temp1 = temp.drop(columns=['Substrate', 'Metabolite','VALUE']).set_index('Product')
    temp1['GENE'] = temp1[['Gene', 'cellname']].apply(lambda x: '_'.join(x), axis=1)
    temp1 = temp1.drop(columns = ['Gene', 'cellname'])
    temp1 = temp1.pivot_table(index = temp1.index, \
                                      columns = 'GENE', \
                                      values = 'log10(Km)')
    temp1 = temp1.T

    # temp 2 contains [metabolites]
    temp2 = temp.drop(columns=['Substrate', 'Metabolite','log10(Km)']).set_index('Product')
    temp2['GENE'] = temp2[['Gene', 'cellname']].apply(lambda x: '_'.join(x), axis=1)
    temp2 = temp2.drop(columns = ['Gene', 'cellname'])
    temp2 = temp2.pivot_table(index = temp2.index, \
                                      columns = 'GENE', \
                                      values = 'VALUE')
    temp2 = temp2.T

    nci_60_final = pd.merge(temp1, temp2, how='inner', left_index=True, right_index=True, suffixes=('_log10(Km)', '_conc'))
    nci_60_final = nci_60_final/100

    # Add everything else:
    # make three human matrix that contains only necessary information
    clean_human_km_average = pd.DataFrame()
    clean_human_km_average['Gene'] = Sabio_RK['Gene names']
    clean_human_km_average['Metabolite'] = Sabio_RK['parameter.associatedSpecies']
    clean_human_km_average['log10(Km)'] = Sabio_RK['log(Km)']
    clean_human_km_average = clean_human_km_average.dropna()

    clean_human_km_median = pd.DataFrame()
    clean_human_km_median['Gene'] = Sabio_RK['Gene names']
    clean_human_km_median['Metabolite'] = Sabio_RK['parameter.associatedSpecies']
    clean_human_km_median['log10(Km)'] = Sabio_RK['log(Km)']
    clean_human_km_median = clean_human_km_median.dropna()

    clean_human_km_binary = pd.DataFrame()
    clean_human_km_binary['Gene'] = Sabio_RK['Gene names']
    clean_human_km_binary['Metabolite'] = Sabio_RK['parameter.associatedSpecies']
    clean_human_km_binary['log10(Km)'] = Sabio_RK['log(Km)']
    clean_human_km_binary = clean_human_km_binary.dropna()

    # make the columns
    metabolites_average = pd.DataFrame(columns=['NAD_average', 'NADH_average', 'NADP_average',\
                                        'NADPH_average', 'GTP_average', 'GDP_average',\
                                        'TTP_average', 'TDP_average', 'CTP_average',\
                                        'CDP_average', 'ATP_average', 'ADP_average',\
                                        'UTP_average', 'UDP_average', 'Acetyl-CoA_average']).fillna(0)
    metabolites_median = pd.DataFrame(columns=['NAD_median', 'NADH_median', 'NADP_median',
                                        'NADPH_median', 'GTP_median', 'GDP_median',
                                        'TTP_median', 'TDP_median', 'CTP_median',
                                        'CDP_median', 'ATP_median', 'ADP_median',
                                        'UTP_median', 'UDP_median', 'Acetyl-CoA_median']).fillna(0)
    metabolites_binary = pd.DataFrame(columns=['NAD_binary', 'NADH_binary', 'NADP_binary',
                                        'NADPH_binary', 'GTP_binary', 'GDP_binary',
                                        'TTP_binary', 'TDP_binary', 'CTP_binary',
                                        'CDP_binary', 'ATP_binary', 'ADP_binary',
                                        'UTP_binary', 'UDP_binary', 'Acetyl-CoA_binary']).fillna(0)

    # construct full matrix
    clean_human_km_average = pd.concat([clean_human_km_average, metabolites_average])
    clean_human_km_median = pd.concat([clean_human_km_median, metabolites_median])
    clean_human_km_binary = pd.concat([clean_human_km_binary, metabolites_binary])

    # binary mapping
    clean_human_km_binary['NAD'] = [1 if ele =='NAD+' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['NADH'] = [1 if ele =='NADH' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['NADP'] = [1 if ele =='NADP+' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['NADPH'] = [1 if ele =='NADPH' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['GTP'] = [1 if ele =='GTP' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['GDP'] = [1 if ele =='GDP' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['TTP'] = [1 if ele =='TTP' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['TDP'] = [1 if ele =='TDP' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['CTP'] = [1 if ele =='CTP' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['CDP'] = [1 if ele =='CDP' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['ATP'] = [1 if ele =='ATP' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['ADP'] = [1 if ele =='ADP' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['UTP'] = [1 if ele =='UTP' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['UDP'] = [1 if ele =='UDP' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary['Acetyl-CoA'] = [1 if ele =='Acetyl-CoA' else 0 for ele in clean_human_km_binary['Metabolite']]
    clean_human_km_binary = clean_human_km_binary.drop(columns = ['Metabolite', 'log10(Km)'])
    clean_human_km_binary = clean_human_km_binary.drop_duplicates(keep='first')
    clean_human_km_binary = clean_human_km_binary.groupby(clean_human_km_binary.Gene).sum()

    # average mapping
    NAD = (clean_human_km_average['Metabolite'] == 'NAD+')
    clean_human_km_average['NAD_average'][NAD] = clean_human_km_average['log10(Km)'][NAD]
    NADH = (clean_human_km_average['Metabolite'] == 'NADH')
    clean_human_km_average['NADH_average'][NADH] = clean_human_km_average['log10(Km)'][NADH]
    NADP = (clean_human_km_average['Metabolite'] == 'NADP+')
    clean_human_km_average['NADP_average'][NADP] = clean_human_km_average['log10(Km)'][NADP]
    NADPH = (clean_human_km_average['Metabolite'] == 'NADPH')
    clean_human_km_average['NADPH_average'][NADPH] = clean_human_km_average['log10(Km)'][NADPH]
    GTP = (clean_human_km_average['Metabolite'] == 'GTP')
    clean_human_km_average['GTP_average'][GTP] = clean_human_km_average['log10(Km)'][GTP]
    GDP = (clean_human_km_average['Metabolite'] == 'GDP')
    clean_human_km_average['GDP_average'][GDP] = clean_human_km_average['log10(Km)'][GDP]
    TTP = (clean_human_km_average['Metabolite'] == 'TTP')
    clean_human_km_average['TTP_average'][TTP] = clean_human_km_average['log10(Km)'][TTP]
    TDP = (clean_human_km_average['Metabolite'] == 'TDP')
    clean_human_km_average['TDP_average'][TDP] = clean_human_km_average['log10(Km)'][TDP]
    CTP = (clean_human_km_average['Metabolite'] == 'CTP')
    clean_human_km_average['CTP_average'][CTP] = clean_human_km_average['log10(Km)'][CTP]
    CDP = (clean_human_km_average['Metabolite'] == 'CDP')
    clean_human_km_average['CDP_average'][CDP] = clean_human_km_average['log10(Km)'][CDP]
    ATP = (clean_human_km_average['Metabolite'] == 'ATP')
    clean_human_km_average['ATP_average'][ATP] = clean_human_km_average['log10(Km)'][ATP]
    ADP = (clean_human_km_average['Metabolite'] == 'ADP')
    clean_human_km_average['ADP_average'][ADP] = clean_human_km_average['log10(Km)'][ADP]
    UTP = (clean_human_km_average['Metabolite'] == 'UTP')
    clean_human_km_average['UTP_average'][UTP] = clean_human_km_average['log10(Km)'][UTP]
    UDP = (clean_human_km_average['Metabolite'] == 'UDP')
    clean_human_km_average['UDP_average'][UDP] = clean_human_km_average['log10(Km)'][UDP]
    Acetyl_CoA = (clean_human_km_average['Metabolite'] == 'Acetyl-CoA')
    clean_human_km_average['Acetyl-CoA_average'][Acetyl_CoA] = clean_human_km_average['log10(Km)'][Acetyl_CoA]
    clean_human_km_average = clean_human_km_average.drop_duplicates(keep='first')
    clean_human_km_average = clean_human_km_average.fillna(value=0)
    clean_human_km_average = clean_human_km_average.groupby(clean_human_km_average.Gene).mean()

    # median mapping
    NAD = (clean_human_km_median['Metabolite'] == 'NAD+')
    clean_human_km_median['NAD_median'][NAD] = clean_human_km_median['log10(Km)'][NAD]
    NADH = (clean_human_km_median['Metabolite'] == 'NADH')
    clean_human_km_median['NADH_median'][NADH] = clean_human_km_median['log10(Km)'][NADH]
    NADP = (clean_human_km_median['Metabolite'] == 'NADP+')
    clean_human_km_median['NADP_median'][NADP] = clean_human_km_median['log10(Km)'][NADP]
    NADPH = (clean_human_km_median['Metabolite'] == 'NADPH')
    clean_human_km_median['NADPH_median'][NADPH] = clean_human_km_median['log10(Km)'][NADPH]
    GTP = (clean_human_km_median['Metabolite'] == 'GTP')
    clean_human_km_median['GTP_median'][GTP] = clean_human_km_median['log10(Km)'][GTP]
    GDP = (clean_human_km_median['Metabolite'] == 'GDP')
    clean_human_km_median['GDP_median'][GDP] = clean_human_km_median['log10(Km)'][GDP]
    TTP = (clean_human_km_median['Metabolite'] == 'TTP')
    clean_human_km_median['TTP_median'][TTP] = clean_human_km_median['log10(Km)'][TTP]
    TDP = (clean_human_km_median['Metabolite'] == 'TDP')
    clean_human_km_median['TDP_median'][TDP] = clean_human_km_median['log10(Km)'][TDP]
    CTP = (clean_human_km_median['Metabolite'] == 'CTP')
    clean_human_km_median['CTP_median'][CTP] = clean_human_km_median['log10(Km)'][CTP]
    CDP = (clean_human_km_median['Metabolite'] == 'CDP')
    clean_human_km_median['CDP_median'][CDP] = clean_human_km_median['log10(Km)'][CDP]
    ATP = (clean_human_km_median['Metabolite'] == 'ATP')
    clean_human_km_median['ATP_median'][ATP] = clean_human_km_median['log10(Km)'][ATP]
    ADP = (clean_human_km_median['Metabolite'] == 'ADP')
    clean_human_km_median['ADP_median'][ADP] = clean_human_km_median['log10(Km)'][ADP]
    UTP = (clean_human_km_median['Metabolite'] == 'UTP')
    clean_human_km_median['UTP_median'][UTP] = clean_human_km_median['log10(Km)'][UTP]
    UDP = (clean_human_km_median['Metabolite'] == 'UDP')
    clean_human_km_median['UDP_median'][UDP] = clean_human_km_median['log10(Km)'][UDP]
    Acetyl_CoA = (clean_human_km_median['Metabolite'] == 'Acetyl-CoA')
    clean_human_km_median['Acetyl-CoA_median'][Acetyl_CoA] = clean_human_km_median['log10(Km)'][Acetyl_CoA]
    clean_human_km_median = clean_human_km_median.drop_duplicates(keep='first')
    clean_human_km_median = clean_human_km_median.fillna(value=0)
    clean_human_km_median = clean_human_km_median.groupby(clean_human_km_median.Gene).median()

    # merge datasets and clean data up
    clean_human_km = pd.merge(clean_human_km_average, clean_human_km_median, how='inner', left_index=True, right_index=True)
    clean_human_km = pd.merge(clean_human_km, clean_human_km_binary, how='inner', left_index=True, right_index=True)
    clean_human_km = clean_human_km.drop(columns = ['log10(Km)_x','log10(Km)_y'])

    os.chdir(input_path)
    for file in os.listdir(os.getcwd()):
        if file.endswith('.csv'):
            filename = os.path.splitext(file)[0]
            train = pd.read_csv(file)
            train['name'] = train['GENE'].str.split('_').str[0]
            train.set_index('name', drop=True, inplace=True)
            train = pd.merge(train, clean_human_km, how='left', left_index=True, right_index=True)
            train.drop(columns=['TCGA_val', 'CNV_val', 'TCGA_annot', 'CNV'], axis=1, inplace=True)
            train = train.replace(0, NaN)
            train = train.dropna(axis=0, how='all', thresh=1)
            train = train.fillna(value=0, axis=0)
            train = train.reset_index().drop(columns='index')
            train = train.set_index('GENE')
            train = pd.merge(train, nci_60_final, how='left', left_index=True, right_index=True)
            train = train.fillna(0)

            csvwriter = filename+'_NCI60.csv'
            train.to_csv(output_path+"./"+csvwriter, encoding='utf-8', index=False)

            return train

input_path = r'C:\Users\scott\Google Drive\Work\UM\Research\Sriram\Projects\Cancer_Modeling\MetOncoFit\raw_data\Krishna'
output_path = r'C:\Users\scott\Google Drive\Work\UM\Research\Sriram\Projects\Cancer_Modeling\MetOncoFit\processed_data\nci60_cofactor_avg_med'
metabolizer(input_path, output_path)
