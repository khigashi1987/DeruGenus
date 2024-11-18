import os
import json
from Bio import Entrez
import time

# NCBI APIの設定
Entrez.email = os.environ.get('ENTREZ_ACCESS_EMAIL')

# Genusリストの読み込み
genera_info = {}
with open('bacteria_archaea_genus_list.txt', 'r') as file:
    for line in file:
        genus_name = line.strip()
        if '(' in genus_name:
            continue
        
        genera_info[genus_name] = {}

        if genus_name.startswith('Candidatus'):
            genera_info[genus_name]['search_name'] = genus_name.split(' ')[1]
        else:
            genera_info[genus_name]['search_name'] = genus_name

print('Total genera:', len(genera_info))

# 論文数カウント
no_hit_genera = []
for i, genus_name in enumerate(genera_info.keys()):
    print(f'{i+1}/{len(genera_info)}: {genus_name}')
    search_name = genera_info[genus_name]['search_name']
    handle = Entrez.esearch(db="pubmed", term=f"{search_name}[Title/Abstract]", retmax=0)
    record = Entrez.read(handle)
    handle.close()
    genera_info[genus_name]['HitCount'] = record['Count']
    print('\tHitCount:', record['Count'])
    
    if record['Count'] == '0':
        no_hit_genera.append(genus_name)

    time.sleep(1)

print('No hit genera:', no_hit_genera)
for genus_name in no_hit_genera:
    del genera_info[genus_name]

# save genera_info as JSON
with open('genera_info.json', 'w') as file:
    json.dump(genera_info, file, indent=4)
