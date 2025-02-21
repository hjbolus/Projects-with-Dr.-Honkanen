#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 09:22:41 2025

@author: harrisbolus
"""
import requests, sys, json, re
import pandas as pd

# Variable regions can be indicated with a period '.' if any amino acid is acceptable at a position, or a bracketed series of acceptable ones.
# An unbracketed single letter abbreviation, or a series of them, is interpreted literally (as an invariant region).
# 'L[ST]PI.E' is equivalent to 'L-[ST]-P-I-x-E' in Prosite notation and matches 'LSPI.E' and 'LTPI.E' where '.' is anything but a space.
# If a single subsequence matches multiple motifs, only the first will be matched. Therefore, it is recommended to order such motifs from most
# to least specific, as shown here.
motifs = [
re.compile(r'L[ST]PI.E'),
re.compile(r'[LCVM][ST]PI.E')
re.compile(r'[LVMIFC][ST]P[ILVM].E')
]

uniprots = ['uniprot1', 'uniprot2', '...']

def get_sequence_from_uniprot(uniprot: str) -> str:
    """Given a Uniprot id, returns the protein's sequence from Uniprot as a string."""
    
    attempts = 0
    params = {"fields": ["sequence"]}
    headers = {"accept": "application/json"}
    base_url = f"https://rest.uniprot.org/uniprotkb/{uniprot}"
    
    response = requests.get(base_url, headers=headers, params=params)
    while response.status_code != 200:
        if attempts < 5:
            attempts+=1
            time.sleep(.2)
            response = requests.get(base_url, headers=headers, params=params)
        else:
            return False
    return response.json()['sequence']['value']

if __name__ == "__main__":
    results = {}
    for uniprot in uniprots:
        print(uniprot)
        if (sequence := get_sequence_from_uniprot(uniprot)):
            for motif in motifs:
                for result in motif.finditer(sequence):
                    position = result.start()
                    pattern = motif.pattern.replace('.','x')
                    match = result[0]
                    if uniprot in results:
                        if not position in results[uniprot]:
                            results[uniprot][position] = (pattern, match)
                    else:
                        results[uniprot] = {position: (pattern, match)}
        else:
            print('Failed to get sequence from uniprot')
