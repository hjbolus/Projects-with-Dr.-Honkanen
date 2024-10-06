#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 18:12:07 2024

@author: harrisbolus
"""

import pandas as pd
import csv
import py4cytoscape as py4
import requests

#--------------parameters-----------------
    #This is the file from a WebGestalt pathway enrichment analysis that says "enrichment_results_wg_result..."
webgestalt_file = '/.../enrichment_results_wg_result ... .txt'

    #These are attributes of the excel file with your proteomics data. Set modsites_column = False if you're using proteomics rather than phosphoproteomics
excel_file = '/.../... .xlsx'
excel_sheet = '...'
modsites_column = False
uniprot_column = '...'
abundance_column = '...'

    #This is just a name that will be used to label output
data_name = '...'

    #Keep 'description' in this list. Otherwise, these should be names of pathways that you think are artifacts of the cell line background.
skiplist = ['description']

#------------define string query function and trend finder-------------------

def getfromstring(method, inputs):                                  #method can be 'get_string_ids' or 'network'
    string_api_url = 'https://string-db.org/api'
    output_format = 'tsv'
    params = {
    'identifiers' : '\r'.join(inputs),
    'species' : 9606, # species NCBI identifier
    'limit' : 1, # only one identifier per input protein
    'echo_query' : 1, # see your input identifiers in the output
    'caller_identity' : 'harris_bolus'}
    request_url = '/'.join([string_api_url, output_format, method])
    results = requests.post(request_url, data=params)
    return(results)

def find_trends(proteinlist, dataset, idcol, valuecol, name):
    direction_list=[]
    for i in proteinlist:
        if all(dataset[dataset[idcol] == i][valuecol] > 0):
            direction_list.append([i,1])
        elif all(dataset[dataset[idcol] == i][valuecol] < 0):
            direction_list.append([i,-1])
        else:
            direction_list.append([i,0])
    return(pd.merge(ourdata, pd.DataFrame(direction_list, columns=[idcol,(name+' trends')]), how='left'))

def fibonacci_of(n):
    if n in {0, 1}:  # Base case
        return n
    return fibonacci_of(n - 1) + fibonacci_of(n - 2)

#---------------put it all together-----------------------------------------

if __name__ == "__main__":
    ourdata = pd.read_excel(excel_file, sheet_name=excel_sheet)
    big_protein_list = set(ourdata[uniprot_column])
    
    with open(webgestalt_file) as csvfile:
        reader = list(csv.reader(csvfile, delimiter='	'))
        prelim_protein_list = [i[-1] for i in reader][1:]
        pathway_list = [i[1] for i in reader]
    
    full_protein_list = []
    for i in prelim_protein_list:
        for j in i.split(';'):
            if j not in full_protein_list:
                full_protein_list.append(j)
    
    mapped_df = pd.DataFrame([{'gene name':i.split('\t')[5], 'uniprotID':i.split('\t')[0], 'id':i.split('\t')[2]} for i in getfromstring('get_string_ids', full_protein_list).text.strip().split("\n")[1:]])
    all_edges = pd.DataFrame([{'source':i.split('\t')[0], 'target':i.split('\t')[1], 'interaction':'interacts'} for i in getfromstring('network', full_protein_list).text.strip().split("\n")[1:]]).drop_duplicates()
    final_nodes = pd.DataFrame()
    
    #skip the first item in the WebGestalt file
    counter = 1
    used_pathway_list = []
    for i in reader:
        if counter > 40:
            break
        if i[1] in skiplist:
            pathway_list.remove(i[1])
            continue
        used_pathway_list.append(i[1])
        counter+=1
            #Get a list of proteins & pick them out from the big dataframe
        protein_list = i[-1].split(';')
        protein_list_df = pd.DataFrame.from_dict({'uniprotID':protein_list, 'pathway':i[1], 'uniqueID':[i[1]+protein for protein in protein_list]})
        mapped_proteins = pd.DataFrame.merge(protein_list_df, mapped_df, on='uniprotID')
    
            #Merge the new dataframe with our data to make the node table
        ourdata = find_trends(full_protein_list, ourdata, uniprot_column, abundance_column, data_name)
        node_table = pd.DataFrame.merge(mapped_proteins, ourdata, how='left', left_on='uniprotID', right_on=uniprot_column)
    
            #Combine IDs with modsite column to make a unique identifier for merging
        if modsites_column:
            node_table.insert(0,'phospho uniprot', node_table['uniprotID'] + ': ' + node_table[modsites_column])
            node_table.insert(0,'phospho gene name', node_table['gene name'] + ': ' + node_table[modsites_column])
    
        final_nodes = pd.concat([final_nodes,node_table])
    
            #Now filter for the edges I need
        preliminary_edge_table = all_edges[all_edges['source'].isin(node_table['id'])]
        edge_table = preliminary_edge_table[preliminary_edge_table['target'].isin(node_table['id'])]
    
            #lastly, send node and edge table to py4cytoscape
        py4.create_network_from_data_frames(node_table, edge_table, title=str(node_table.loc[0,'pathway']))#+ ' grouped')
            #Group them if you want individual portraits
        #py4.groups.create_group(group_name=i[1], nodes='all')
    
    #merge using unique IDs to preserve each pathway diagram
    try:
        py4.tools.merge_networks(sources=used_pathway_list, operation='union', node_keys=['uniqueID' for column in used_pathway_list])
    except TypeError:
        pass
    
    # group each pathway
    for i in used_pathway_list:
        py4.create_group_by_column(i, column='pathway', value=i)
    
    #arrange groups
    node_table = py4.commands.cyrest_post(operation='commands/node/get attribute')['data']
    
    for group in py4.list_groups()['groups']:
        group = py4.get_group_info(group)['name']
        print(group)
        protein_nodes = [f'uniqueID:{i["uniqueID"]}' for i in node_table if (i['pathway'] == group and i['uniqueID'])]
        py4.commands.cyrest_post(operation='commands/layout/attribute-circle', body={'nodeList':','.join(protein_nodes)+',gene name:'+group, 'spacing':'0'})
        start = 0
        layer = 1
        while start < len(protein_nodes):
            diff = fibonacci_of(layer+6)
            end = start+diff
            if end > len(protein_nodes):
                end = len(protein_nodes)
            if end == len(protein_nodes)-1:
                end+=1
            n = end-start
            correction = 1+(layer*0.1)
            if diff != n:
                if n == 3:
                    correction = 1.8
                elif n ==4 or n==5:
                    correction = 0.9
                elif n >= 6:
                    correction = 1
                    if n < 9:
                        j = 0
                    elif n < 16:
                        j = 1
                    elif n < 25:
                        j = 3
                    elif n < 36:
                        j = 5
                    for k in [1.5, 1.4, 1.3, 1.2, 1.1, 1.05, 1.02, 1.008, 1.005][j:j+layer]:
                        correction = correction*k
            py4.commands.cyrest_post(operation='commands/layout/attribute-circle', body={'nodeList':','.join(protein_nodes[start:end]), 'spacing':70*correction})
            diff+=1
            layer+=1
            start=end
