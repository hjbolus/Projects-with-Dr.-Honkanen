#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 08:07:39 2024

@author: harrisbolus
"""

import requests
import pandas as pd
import time
from functools import lru_cache

uniprot_list = ['uniprot1', 'uniprot2', '...']

#Variable regions can be indicated with an 'x' if any amino acid is acceptable at a position, or a bracketed series of acceptable ones.
#An unbracketed single letter abbreviation, or a series of them, is interpreted literally (as an invariant region), and should be separated from variable regions by a dash.
motifs = ['LSPI-x-E', 
          '[LCVM]-SPI-x-E', 
          '[LVMIFC]-SP-[ILVM]-x-E',
          'L-TPI-x-E', 
          '[LCVM]-TPI-x-E', 
          '[LCVMIF]-TP-[ILVM]-x-E'
          ]

def test_motif_element(characters: str, motif_element: str) -> bool:
    if motif_element == 'x':
        return True
    elif motif_element[0] == '[':
        assert len(characters)==1 and len(motif_element)>2, print(f'characters {characters} is too long to match the single motif range, or motif element {motif_element} is too short')
        return any(characters==element for element in motif_element[1:-1])
    else:
        try:
            return all(characters[i]==motif_element[i] for i in range(len(motif_element)))
        except:
            print(f' checking if {characters} match {motif_element}')

def fits_motif(frame: str, motif: list) -> bool:
    assert isinstance(motif, list)
    assert isinstance(frame, str)
    if len(frame) == sum([len(i) if i[0] !='[' else 1 for i in motif]):
        matches = []
        position = 0
        for element in motif:
            
            if element[0] == '[':
                element_length = 1
            elif element == 'x':
                position+=1
                matches.append(True)
                continue
            else:
                element_length = len(element)
            matches.append(test_motif_element(frame[position:position+element_length], element))
            position += element_length
        return all(matches)
    else:
        return False        

def search_for_motif_in_full_sequence(sequence: str, motif: list) -> tuple:
    solid_elements = [i for i in motif if i[0] not in {'[', ']', 'x'}]
    if all(element in sequence for element in solid_elements):
        motif_length = sum([len(i) if i[0] !='[' else 1 for i in motif])
        sequence_length = len(sequence)
        for i in range(sequence_length):
            if sequence_length-i>=motif_length:
                frame = sequence[i:i+motif_length]
                if fits_motif(frame, motif):
                    return (frame, i)
    return False

@lru_cache(maxsize=len(motifs))
def get_reference_element(motif: tuple) -> str:
    assert isinstance(motif, tuple)
    longest = ('', 0)
    for i in range(len(motif)):
        element = motif[i]
        if element[0] not in {'[', ']', 'x'}:
            if len(element) > len(longest[0]):
                longest = (element, i)
    return longest

def find_positions_of_reference_element_in_sequence(reference_element_tuple: tuple, sequence: str) -> list:
    reference_element, reference_position_in_motif = reference_element_tuple
    reference_element_len = len(reference_element)
    
    reference_positions_in_sequence = []
    if reference_element and reference_element in sequence:
        pieces = sequence.split(reference_element)
        read_frame = 0
        for piece in pieces[:-1]:
            read_frame+=len(piece)
            reference_positions_in_sequence.append(read_frame+1)
            read_frame+=reference_element_len
        return reference_positions_in_sequence
    else:
        return []

def find_frames(sequence: str, motif: list) -> list:
    reference_element, reference_position_in_motif = get_reference_element(tuple(motif))
    frames = []
    if reference_element[0]:
        candidate_positions = find_positions_of_reference_element_in_sequence((reference_element, reference_position_in_motif), sequence)
        for candidate_position in candidate_positions:
            start = (candidate_position - reference_position_in_motif) - 1
            end = start + sum([len(i) if i[0] !='[' else 1 for i in motif])
            frames.append((sequence[start:end], start))
        return frames
    else:
        return None
    
def search_for_motif(sequence:str, motif:str) -> list:
    motif = motif.split('-')
    frames = find_frames(sequence, motif)
    if frames:
        instances = []
        for frame in frames:
            if fits_motif(frame[0], motif):
                instances.append(frame)
        return instances if instances else False
    else:
        return False

proteins = {}
for uniprot in uniprot_list:
    print(uniprot)
    time.sleep(0.1)
    attempts = 0
    response = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniprot}.txt')

    while response.status_code != 200:
        if attempts < 5:
            attempts+=1
            time.sleep(.2)
            response = requests.get(f'https://rest.uniprot.org/uniprotkb/{uniprot}.txt')
        else:
            continue
    
    sequence = response.text.split(';')[-1].replace('\n','').replace(' ','').replace('/','')
    for motif in motifs:
        results = search_for_motif(sequence, motif)
        if results:
            for result in results:
                position = result[-1]
                if uniprot in proteins:
                    if not position in proteins[uniprot]:
                        proteins[uniprot][position] = (motif, result[0])
                else:
                    proteins[uniprot] = {position: (motif, result[0])}
