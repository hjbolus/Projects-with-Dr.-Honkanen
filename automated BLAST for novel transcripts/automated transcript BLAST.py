#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Performs a BLASTn using the NCBI core_nt database for a list of transcripts.
"""
Created on Sat Oct 19 11:03:45 2024

@author: harrisbolus
"""

import requests
from Bio import Blast
from Bio.Blast import NCBIWWW
import zipfile
from datetime import datetime
from time import sleep
import pandas as pd

#---------------------------Parameters-------------------------------------------------------

# Ensembl transcript IDs
transcript_list = ['ENST00000623324', 'ENST00000648868', ...]

output_location = '/Users/.../.../'

# Would you like to BLAST the sense or antisense sequence of the transcripts above?
antisense = False

Blast.email = 'your email'
NCBIWWW.email = Blast.email

# Dictionary of completed exons from given transcripts. Generated by generate_dictionary_of_completed_transcripts()
done = {}

#---------------------------------------------------------------------------------------------

nucleotide_pairs = {'G':'C','C':'G','A':'T','T':'A'}

is_antisense = ' antisense' if antisense else ''

def get_exon_sequences(transcript_id: str) -> dict:
    # Get exon IDs
    lookup_url = f"https://rest.ensembl.org/lookup/id/{transcript_id}?expand=1"
    headers = {"Content-Type": "application/json"}
    response = requests.get(lookup_url, headers=headers)
    response.raise_for_status()
    transcript_data = response.json()

    exon_ids = [exon['id'] for exon in transcript_data['Exon']]

    # Get exon sequences
    exon_sequences = {}
    for exon_id in exon_ids:
        sequence_url = f"https://rest.ensembl.org/sequence/id/{exon_id}"
        response = requests.get(sequence_url, headers=headers)
        response.raise_for_status()
        exon_data = response.json()
        exon_sequences[exon_id] = exon_data['seq']

    return exon_sequences

def antisense(sequence: str) -> str:
    """returns the reverse complement of the given sequence"""
    return ''.join([nucleotide_pairs[i] for i in sequence[::-1]])

def generate_dictionary_of_completed_transcripts(filenames: str) -> dict:
    """Allows you to skip previously BLASTed transcripts without editing your source file. 
    If you get a timeout or other error, you can copy and paste the names of the files 
    previously generated into a single string separated by \n. Copy the dict and set it 
    as the value of the variable "done" above."""
    a = filenames.replace('.zip','').replace(' exon','')
    b = [i.split(' ') for i in a.split('\n')]
    c = set(i[0] for i in b)
    dict1 = {}
    for i in c:
        dict1[i] = [j[1] for j in b if j[0] == i]
    return dict1

if __name__ == "__main__":
    for transcript_id in transcript_list:
        print(f'Working through exons for transcript {transcript_id}')
        exon_sequences = get_exon_sequences(transcript_id)
        counter = 1
        for ensembl_id, exon_seq in list(exon_sequences.items()):
            if done:
                if transcript_id in done: #get remaining exons for this id
                    if str(counter) in done[transcript_id]:
                        print(f'skipping {transcript_id} exon {counter}')
                        counter+=1
                        continue
            if antisense:
                exon_seq = antisense(exon_seq)
            print(f'requesting exon {counter}')
            start_time = datetime.now()
            result = Blast.qblast('blastn', 'core_nt', exon_seq, format_type='JSON2', entrez_query='Homo sapiens[organism]').read()
            print('got it')
    
            with open(f'{output_location}{transcript_id} exon {counter}{is_antisense}.zip', 'wb') as out_stream:
                out_stream.write(result)
            counter+=1
            time_elapsed = (datetime.now() - start_time).total_seconds()
            if time_elapsed < 15:
                sleep(15-time_elapsed)
