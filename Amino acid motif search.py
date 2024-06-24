import requests
import pandas as pd
import time

motifs = ['LSPI-x-E', '[LCVM]-SPI-x-E', '[LVMIFC]-SP-[ILVM]-x-E', 'L-TPI-x-E', '[LCVM]-TPI-x-E', '[LCVMIF]-TP-[ILVM]-x-E'] #SLiMs
file = '/.../.../... .xlsx'
uniprot_column = ...

df = pd.read_excel(file)
proteins = df[uniprot_column] #column with Uniprot IDs

def test_motif_element(characters, motif_element):
  #boolean test for a match between a given character and a motif element
    if motif_element == 'x':
        return True
    elif motif_element[0] == '[':
        assert len(characters)==1, print(characters, motif_element)
        return any(characters==element for element in motif_element[1:-1])
    else:
        return all(characters[i]==motif_element[i] for i in range(len(motif_element)))

def search_for_motif(sequence, motif):
  #searches for a complete match (according to test_motif_element()) in any frame of a given sequence. 
  #If found, returns the position of the first AA in the motif. This returns only the first occurrence, 
  #but the position can be used to search again on the remainder of the sequence if desired.

    motif = motif.split('-')
    motif_length = sum([len(i) if i[0] !='[' else 1 for i in motif])
    sequence_length = len(sequence)
    for i in range(sequence_length):
        list1 = []
        if sequence_length-i>=motif_length:
            j = i
            for motif_element in motif:
                j+= len(motif_element) if motif_element[0] != '[' else 1
                list1.append(test_motif_element(sequence[i:j], motif_element))
                i=j
            if all(list1):
                return i-motif_length
    return False


found_motifs = {protein: [] for protein in proteins}
for protein in proteins:
    time.sleep(0.1)
    response = requests.get(f'https://rest.uniprot.org/uniprotkb/{protein}.txt')
    if response.status_code == 200:
        sequence = response.text.split('SEQUENCE')[-1].split(';')[-1].replace('\n','').replace(' ','').replace('/','')
        for motif in motifs:
            if position:=search_for_motif(sequence, motif):
                found_motifs[protein].append((position, motif))
                print(protein, position, motif)
                break #remove if a comprehensive search is desired
    else:
        print(response.status_code)

new_df = pd.DataFrame(found_motifs, columns=[uniprot_column,'motif'])
df.merge(new_df, how='left', on=uniprot_column).to_excel('/'.join(file.split('/')[0:-1]) + file.split('/')[-1].split('.')[0]+'motif search.xslx')

