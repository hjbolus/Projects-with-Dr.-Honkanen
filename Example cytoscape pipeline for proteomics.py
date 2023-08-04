import pandas as pd
import csv
import py4cytoscape as py4
import requests

#--------------parameters-----------------
    #This is the file from a WebGestalt pathway enrichment analysis that says "enrichment_results_wg_result..."
webgestalt_file = '/Users/harrisbolus/Desktop/Research/Dr. Honkanen/PPP5c KO lab/PPP5C KO Figure prep/Webgestalt downloads/2 prot run 2/enrichment_results_wg_result1683083193.txt'

    #These are attributes of the excel file with your proteomics data. Set modsites_column = False if you're using proteomics rather than phosphoproteomics
excel_file = '/Users/harrisbolus/Desktop/Research/Dr. Honkanen/PPP5c KO lab/HB 05072023 PPP5C KO proteomics vs mouse comparison.xlsx'
excel_sheet = '2 prot run 2'
modsites_column = False
uniprot_column = 'accession'
abundance_column = 'average'

    #This is just a name that will be used to label output
data_name = 'PPP5C prot run 2'

    #Keep 'description' in this list. Otherwise, these should be names of pathways that you think are artifacts of the cell line background.
skiplist = ['description', 'Human immunodeficiency virus 1 infection', 'Prostate cancer','Human T-cell leukemia virus 1 infection','Pathogenic Escherichia coli infection','Viral carcinogenesis','Alcoholism','Bacterial invasion of epithelial cells','Endometrial cancer','Salmonella infection','Proteoglycans in cancer','Hepatocellular carcinoma','Herpes simplex infection','Choline metabolism in cancer']

#------------define string query function and trend finder-------------------

def getfromstring(method, inputs):                                  #method can be 'get_string_ids' or 'network'
    string_api_url = 'https://version-11-5.string-db.org/api'
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

#---------------put it all together-----------------------------------------

reader = csv.reader(open(webgestalt_file),delimiter='	')
pathway_list = [i[1] for i in reader]

reader = csv.reader(open(webgestalt_file),delimiter='	')
prelim_protein_list = [i[-1] for i in reader][1:]

ourdata = pd.read_excel(excel_file, sheet_name=excel_sheet)
big_protein_list = set(ourdata[uniprot_column])

real_protein_list = []
for i in prelim_protein_list:
    for j in i.split(';'):
        if j not in real_protein_list:
            real_protein_list.append(j)

big_mapped_df = pd.DataFrame([{'gene name':i.split('\t')[5], 'uniprotID':i.split('\t')[0], 'id':i.split('\t')[2]} for i in getfromstring('get_string_ids', real_protein_list).text.strip().split("\n")[1:]])
all_edges = pd.DataFrame([{'source':i.split('\t')[0], 'target':i.split('\t')[1], 'interaction':'interacts'} for i in getfromstring('network', real_protein_list).text.strip().split("\n")[1:]]).drop_duplicates()
final_nodes = pd.DataFrame()

reader = csv.reader(open(webgestalt_file),delimiter='	')
#skip the first item in the WebGestalt file
for i in reader:
    if i[1] in skiplist:
        pathway_list.remove(i[1])
        continue

        #Get a list of proteins & pick them out from the big dataframe
    protein_list = i[-1].split(';')
    protein_list_df = pd.DataFrame.from_dict({'uniprotID':protein_list, 'pathway':i[1], 'uniqueID':[i[1]+protein for protein in protein_list]})
    mapped_proteins = pd.DataFrame.merge(protein_list_df, big_mapped_df, on='uniprotID')

        #Merge the new dataframe with our data to make the node table
    ourdata = find_trends(real_protein_list, ourdata, uniprot_column, abundance_column, data_name)
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
    py4.tools.merge_networks(sources=pathway_list, operation='union', node_keys=['uniqueID' for column in pathway_list])
except TypeError:
    pass

#group each pathway
for i in pathway_list:
    py4.create_group_by_column(i, column='pathway', value=i)

#Add edges between pathways
group_list = pd.DataFrame([py4.get_group_info(i) for i in py4.list_groups()['groups']])

pathway_edges = pd.merge(final_nodes.loc[:,['pathway','id']], all_edges, how='inner', left_on='id', right_on='source').drop('id', axis=1)
pathway_edges = pd.merge(final_nodes.loc[:,['pathway','id']], pathway_edges, how='inner', left_on='id', right_on='target', suffixes=('_target', '_source')).drop_duplicates(subset=['pathway_source','pathway_target']).drop(['id','source','target'], axis=1)
pathway_edges = pathway_edges[pathway_edges['pathway_source']!=pathway_edges['pathway_target']]

done = set()
for i in pathway_edges.values:
    if (i[0],i[1]) not in done:
        py4.networks.add_cy_edges([str(group_list[group_list['name'] == i[0]]['group'].values[0]),str(group_list[group_list['name'] == i[1]]['group'].values[0])])
        done.add((i[1],i[0]))
