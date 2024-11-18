import os
import requests
import json
from Bio import Entrez
import time

# Base URL for Europe PMC API
base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# NCBI APIの設定
Entrez.email = os.environ.get('ENTREZ_ACCESS_EMAIL')

# Load genera_info
with open('genera_info.json', 'r') as file:
    genera_info = json.load(file)

# For top 1000 genera by HitCount
top_genera = sorted(genera_info.items(), key=lambda x: int(x[1]['HitCount']), reverse=True)[:1000]

# 各Genusのトップ20論文を取得
for i, g in enumerate(top_genera):

    if os.path.exists(f'./abstracts/abstracts_{i}.json'):
        print(f'{i+1}/{len(top_genera)}: {g[0]} (already exists)')
        continue

    genus_name = g[0]
    search_name = g[1]['search_name']
    hit_count = g[1]['HitCount']

    genus_info = {
        'Genus': genus_name,
        'PubMed Hit Count': hit_count,
        'Papers': []
    }

    print(f'{i+1}/{len(top_genera)}: {genus_name}')

    # search articles with top cited
    response = requests.get(f"{base_url}?query=ABSTRACT:{search_name}%20sort_cited:y&format=json&pageSize=20")
    if response.status_code == 200:
        data = response.json()
        for result in data.get("resultList", {}).get("result", []):
            title = result.get("title", "No Title")
            citation_count = result.get("citedByCount", 0)
            pmid = result.get("pmid", "N/A")
            print(f"\tPMID: {pmid}, Title: {title}, Citations: {citation_count}")
            genus_info['Papers'].append({'PMID': pmid, 'Title': title, 'Citation Count': citation_count})
    else:
        print(f"\tError: {response.status_code}")

    # Abstractを取得
    for pind in range(len(genus_info['Papers'])):
        pmid = genus_info['Papers'][pind]['PMID']
        if not pmid.isdigit():
            genus_info['Papers'][pind]['Abstract'] = ""
            continue

        try:
            fetch_handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            fetch_record = Entrez.read(fetch_handle)
            fetch_handle.close()
            abstract = fetch_record['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
            abstract = " ".join(abstract)
            time.sleep(1)
        except KeyError:
            abstract = ""

        genus_info['Papers'][pind]['Abstract'] = abstract
    
    # Save as JSON
    with open(f'./abstracts/abstracts_{i}.json', 'w') as file:
        json.dump(genus_info, file, indent=4)