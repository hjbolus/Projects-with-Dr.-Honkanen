#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 12:47:16 2024

@author: harrisbolus
"""
    # statistics of BLAST scores: https://www.ncbi.nlm.nih.gov/BLAST/tutorial/

import zipfile
import os
import json
import pandas as pd

# List of directories to parse BLAST files from
directories = ['/Users/.../.../']

# What do you want to name the excel file?
output_file = '/Users/.../... .xlsx'

def pull_names_out_from_parentheses(string: str) -> list:
    """Given a string containing one or more sets of matched single-layer parentheses, returns a list of their contents"""
    split_string = string.split('(')[1:]
    return set([i.split(')')[0] for i in split_string])
    
row_list = []
for directory in directories:
    print(directory)
    os.chdir(directory)
    zip_file_list = os.listdir()

    for file_name in zip_file_list:
        if file_name == '.DS_Store':
            continue
        transcript, suffix = file_name.split(' exon ')
        exon_number, suffix = suffix.split('.')

        zip_file = zipfile.ZipFile(file_name)
        results_json = [i for i in zip_file.namelist() if '_1' in i]
        assert len(results_json) == 1
        stream = zip_file.open(results_json[0])
        data = json.loads(stream.read().decode())
        hit_list = data['BlastOutput2']['report']['results']['search']['hits']

        for entry in hit_list:
            descriptions = [i for i in entry['description'] if \
                                            ('Homo sapiens isolate ' not in i['title'] and \
                                            'clone' not in i['title'] and \
                                            'ynthetic' not in i['title'] and \
                                            'Homo sapiens uncharacterized lncRNA' not in i['title'] and \
                                            'human bac library' not in i['title'].lower() and \
                                            'human pac library' not in i['title'].lower() and \
                                            'lncRNA gene, complete sequence' not in i['title'] and \
                                            'lncRNA for PREDICTED lncRNA' not in i['title'] and \
                                            'long intergenic non-protein coding RNA' not in i['title'] and \
                                            'Homo sapiens uncharacterized protein (LOC' not in i['title'])
                            ]
            if not descriptions:
                continue
            hsps = entry['hsps']

            plus_strand_evalues = [i['evalue'] for i in hsps if i['hit_strand'] == 'Plus' and i['evalue'] < 0.05]
            if plus_strand_evalues:
                top_plus_strand_hit = [i for i in hsps if i['hit_strand'] == 'Plus' and i['evalue'] == min(plus_strand_evalues)][0]
            else:
                top_plus_strand_hit = False

            minus_strand_evalues = [i['evalue'] for i in hsps if i['hit_strand'] == 'Minus' and i['evalue'] < 0.05]
            if minus_strand_evalues:
                top_minus_strand_hit = [i for i in hsps if i['hit_strand'] == 'Minus' and i['evalue'] == min(minus_strand_evalues)][0]
            else:
                top_minus_strand_hit = False
                
            if top_plus_strand_hit or top_minus_strand_hit:
                row = {'transcript': transcript,
                   'exon number': exon_number,
                   'index': entry['num'],
                   
                   # descriptions
                   'refseq id': ', '.join([i['accession'] for i in descriptions]),
                   'description': (description:=' | '.join([i['title'] for i in descriptions if (' ' not in i['title'] or 'human' not in i['title'])])),
                   'gene names': ', '.join(pull_names_out_from_parentheses(description)),

                   # top plus strand hit
                   'bit_score (plus strand)': top_plus_strand_hit['bit_score'] if top_plus_strand_hit else 'NA',
                   'score (plus strand)': top_plus_strand_hit['score'] if top_plus_strand_hit else 'NA',
                   'e-value (plus strand)': top_plus_strand_hit['evalue'] if top_plus_strand_hit else 'NA',
                   'gaps (plus strand)': top_plus_strand_hit['gaps'] if top_plus_strand_hit else 'NA',
                   'align length (plus strand)': top_plus_strand_hit['align_len'] if top_plus_strand_hit else 'NA',

                   # top minus strand hit
                   'bit_score (minus strand)': top_minus_strand_hit['bit_score'] if top_minus_strand_hit else 'NA',
                   'score (minus strand)': top_minus_strand_hit['score'] if top_minus_strand_hit else 'NA',
                   'e-value (minus strand)': top_minus_strand_hit['evalue'] if top_minus_strand_hit else 'NA',
                   'gaps (minus strand)': top_minus_strand_hit['gaps'] if top_minus_strand_hit else 'NA',
                   'align length (minus strand)': top_minus_strand_hit['align_len'] if top_minus_strand_hit else 'NA'
                   }

                row_list.append(row)
df = pd.DataFrame(row_list)
df.to_excel(output_file)
