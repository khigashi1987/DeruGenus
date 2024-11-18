import os
import pandas as pd

DB_DIR = '/Users/higashi/Desktop/DB/NCBITaxonomy'

nodes_file = os.path.join(DB_DIR, 'nodes.dmp')
names_file = os.path.join(DB_DIR, 'names.dmp')

nodes = pd.read_csv(nodes_file, sep='|', header=None, usecols=[0, 1, 2], 
                    names=["tax_id", "parent_tax_id", "rank"], 
                    skipinitialspace=True)
nodes['tax_id'] = nodes['tax_id'].astype(int)
nodes['parent_tax_id'] = nodes['parent_tax_id'].astype(int)
# 前後のタブ文字を削除
nodes = nodes.applymap(lambda x: x.strip() if isinstance(x, str) else x)

names = pd.read_csv(names_file, sep='|', header=None, usecols=[0, 1, 3], 
                    names=["tax_id", "name_txt", "name_class"], 
                    skipinitialspace=True)
# 前後のタブ文字を削除
names = names.applymap(lambda x: x.strip() if isinstance(x, str) else x)
# 学術名のみ
names = names[names['name_class'] == 'scientific name']

# BacteriaとArchaeaのtax_id
bacteria_tax_id = 2
archaea_tax_id = 2157

# Bacteria/Archaeaの子孫ノードを取得する関数
def get_descendants(tax_id, nodes_df):
    descendants = set()
    stack = [tax_id]
    while stack:
        current = stack.pop()
        descendants.add(current)
        children = nodes_df[nodes_df['parent_tax_id'] == current]['tax_id'].tolist()
        stack.extend(children)
    return descendants

# Bacteria/Archaeaの子孫ノードを取得
bacteria_descendants = get_descendants(bacteria_tax_id, nodes)
archaea_descendants = get_descendants(archaea_tax_id, nodes)

# Bacteria/Archaeaの全ノード
target_tax_ids = bacteria_descendants.union(archaea_descendants)

# Genusに限定
genus_nodes = nodes[(nodes['tax_id'].isin(target_tax_ids)) & (nodes['rank'] == 'genus')]

# Genus名を取得
genus_names = names[names['tax_id'].isin(genus_nodes['tax_id'])]

# 結果を保存
genus_names['name_txt'].to_csv('bacteria_archaea_genus_list.txt', index=False, header=False)