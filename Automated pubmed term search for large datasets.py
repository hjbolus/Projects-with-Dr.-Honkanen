import requests
import time
import pandas as pd
from datetime import datetime
file = '/.../.../... .xlsx'
df = pd.read_excel(file).drop_duplicates()

base_url = 'https://pubmed.ncbi.nlm.nih.gov/?term='

if __name__ == "__main__":
    df1 = df.to_dict('index')
    # For completing partial files
    # df1 = {i:df1[i] for i in df1 if type(df1[i]['pubmed results']) == float}
    
    length = len(df)
    for i in df1.keys():
        print(f'{i}/{length}')
        row = df1[i]
        modsites = row['modsites'].split(':')
        for modsite in modsites:
            query = base_url+row['Gene name']+'+'+modsite
            time.sleep(0.05)
            try:
                response = requests.get(query)
            except TimeoutError:
                time.sleep(30)
                try:
                    response = requests.get(query)
                except TimeoutError:
                    result = 'try again'
    
            if response.status_code == 200:
                if  'The following term was not found in PubMed' in response.text or 'Your search was processed without automatic term mapping because it retrieved zero results' in response.text:
                    print('False', query)
                    result = 'No results'
                else:
                    print('True', query)
                    result = query
    
                if result:
                    if 'pubmed results' in row.keys() and not type(row['pubmed results']) == float:
                        row['pubmed results'] = row['pubmed results']+'; ' + result
                    else:
                        row['pubmed results'] = result
            else:
                print(response.status_code, '\n')
    today = datetime.now()
    pd.DataFrame.from_dict(df1, orient='index').to_excel(file.split('.')[0] + f'{today[1]+today[2]+today[0]}.xlsx')
