import pandas as pd
from collections import *
from tqdm import tqdm
from ete3 import Tree
from glob import glob

g2file = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1783paths.tsv',sep='\t',index_col=0)

target_kos = [
"K04751",
"K04752",
"K13599",
"K13598",
"K07708",
"K03320",
"K00990",]

g2pop = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1783Ruegeria_MCs.tsv',sep='\t',index_col=0)
g2pop = g2pop['MC'].to_dict()
pop2g = defaultdict(list)
for g,pop in g2pop.items():
    pop2g[pop].append(g)
    
kegg_df = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/annotations/KEGG_anno.tab',sep='\t',index_col=0,low_memory=False)
kegg_df2 = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/data_processing/annotations/KEGG_latest/kegg_anno.tsv',sep='\t',index_col=0)
for gid in tqdm(g2file.index[g2file['nanopore']=='YES']):
    if gid not in kegg_df.index or gid not in kegg_df2.index:
        continue
    kegg_df.loc[gid,:] = pd.np.nan
    kegg_df.loc[gid,kegg_df2.columns]=kegg_df2.loc[gid,:]
kegg_df = kegg_df.loc[[_ for _ in g2pop if _ in kegg_df.index],:]

num_df = kegg_df.fillna('').applymap(lambda x: 0 if x=='' else len(x.split(',')))

genomes1 = [g for pop in ['MC6','MC225','MC223','MC13','MC156','MC87'] for g in pop2g.get(pop,[]) ]
genomes2 = set(num_df.index).difference(genomes1)
genomes1 = [_ for _ in genomes1 if _ in num_df.index]
genomes2 = [_ for _ in genomes2 if _ in num_df.index]
print(len(genomes1),len(genomes2))

g2metadata = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1837Ruegeria_MCs.tsv',sep='\t',index_col=0)
coral_gs = list(g2metadata.index[g2metadata['Ecosystem']=='coral'])
coral_gs1 = [_ for _ in genomes1 if _ in coral_gs]


def get_copynumber(tre):
    leaf = tre.get_leaf_names()
    genome2num = defaultdict(int)
    for l in leaf:
        genome2num[l.split('_')[0]]+=1
    return round(np.mean(list(genome2num.values())),2)

ko2tre = {}
tpath2n = {}
ko2inter = {}
for tpath in tqdm(glob('/mnt/ivy/thliao/project/coral_ruegeria/data_processing/denitrification_papaers/KO_phylogeny/*.newick')):
    ko = tpath.split('/')[-1].replace('.newick','')
    if 'coral' in ko:
        continue
    try:
        tre = Tree(tpath)
    except:
        continue
    tre.set_outgroup(tre.get_midpoint_outgroup())
    ko2tre[ko] = tre
    largest_ratio = 0
    for n in tre.traverse():
        all_leafs = n.get_leaf_names()
        genomes = set([_.split('_')[0] for _ in all_leafs])
        genomes = [_ for _ in genomes if _ in coral_gs]
        if len(genomes)==0:
            continue
        intersec = [_ for _ in genomes if _ in genomes1]
        #ra1,ra2 = len(intersec)/len(genomes1),
        if len(intersec)/len(genomes)>=0.7:
            ratio = len(intersec)/len(coral_gs1)
            if ratio >=largest_ratio:
                tpath2n[ko] = (n,ratio)
                largest_ratio = ratio
                ko2inter[ko] = [_ for _ in all_leafs if _.split('_')[0] in intersec]


ko_info = pd.read_csv('/home-user/thliao/db/protein_db/kegg/v20230301/ko_list', sep='\t', index_col=0)
ko2info = {ko: i for ko, i in ko_info['definition'].to_dict().items()}

target_kos = []
for ko,v in sorted(tpath2n.items(),key=lambda x:x[1][1]):
    if v[1] >= 0.8 and get_copynumber(ko2tre[ko])<=2.3:
        print(ko2info[ko],ko,v[1])
        target_kos.append(ko)


n = tre.get_common_ancestor('AE2,L1'.split(','))

tre = Tree(f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_Speciestree.newick',format=3)
gids = list(tre.get_leaf_names())
gids.remove('FO2')
from bin.format_newick import renamed_tree

tre = renamed_tree(tre)
tre.write(outfile='/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_Speciestree_named.newick',format=3)

sids = list(tre.get_leaf_names())

sub_pop2g = defaultdict(list)
for g,pop in g2pop.items():
    if g in tre.get_leaf_names():
        sub_pop2g[pop].append(g)

def get_mc_ratio(leaf):
    mc_ = set([g2pop[_] for _ in leaf])
    for m in mc_:
        full = sub_pop2g[m]
        if len(full)<=3:
            continue
        intersec = set(leaf).intersection(full)
        #diff = set(leaf).intersection(full)
        #print(m,len(intersec)/len(full))
        if len(intersec)/len(full)>0.9 and len(intersec)/len(leaf)>0.8:
            return (m,len(leaf))
    return None,None


mc2size = {}
mc2node = {}
for n in tre.traverse():
    mc,size = get_mc_ratio(n.get_leaf_names())
    if mc is None:
        continue
    if mc not in mc2node:
        mc2node[mc]=n
        mc2size[mc] = size
    else:
        if size >= mc2size.get(mc,0):
            mc2node[mc] = n
            mc2size[mc] = size


with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_Speciestree_named.collasped.txt','w') as f1:
    f1.write('COLLAPSE\nDATA\n')
    for m,v in mc2node.items():
        f1.write(f"{v.name}\n")


with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_Speciestree_named.labels.txt','w') as f1:
    f1.write('LABELS\nSEPARATOR TAB\nDATA\n')
    for m,v in mc2node.items():
        f1.write(f"{v.name}\t{m}\n")

ko2mc2node = defaultdict(dict)
for ko,kotre in tqdm(ko2tre.items()):
    if tpath2n[ko][1] < 0.8:
        continue
    _mc2size = {}
    _mc2node = {}    
    for n in kotre.traverse():
        leaf = set([_.split('_')[0] for _ in n.get_leaf_names() if _.split('_')[0] in sids])
        mc,size = get_mc_ratio(leaf)
        #print(mc,size)
        if mc is None:
            continue
        if mc not in _mc2node:
            _mc2node[mc]=n
            _mc2size[mc] = size
        else:
            if size >= _mc2size.get(mc,0):
                _mc2node[mc] = n
                _mc2size[mc] = size
    ko2mc2node[ko] = _mc2node


ko = 'K03320'

mc1='MC6'
mc2='MC13'
tre = ko2tre[ko]
copy_tre = tre.copy()
for i in copy_tre.get_leaves():
    i.name = i.name.split('_')[0]
copy_tre.write(outfile=f'./tmp_{ko}.newick')

n1,n2 = ko2mc2node[ko][mc1],ko2mc2node[ko][mc2]
mn = tpath2n[ko][0]
print(n1 in list(mn.traverse()),n2 in list(mn.traverse()))

tre_dis = tre.get_distance(ko2mc2node[ko][mc1],ko2mc2node[ko][mc2])/tre.get_farthest_leaf()[1]*100
spe_tre = Tree('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_Speciestree_named.newick',format=3)
spe_tredis = spe_tre.get_distance(mc2node[mc1],mc2node[mc2])/spe_tre.get_farthest_leaf()[1]*100
print(tre_dis,spe_tredis)


for ko in target_kos:
    
    mc1='MC6'
    mc2='MC13'
    mn = tpath2n[ko][0]
    if mc2 in ko2mc2node[ko] and mc1 in ko2mc2node[ko]:
        n1,n2 = ko2mc2node[ko][mc1],ko2mc2node[ko][mc2]
        tre = ko2tre[ko]
        tre_dis = tre.get_distance(ko2mc2node[ko][mc1],ko2mc2node[ko][mc2])/tre.get_farthest_leaf()[1]*100
        spe_tre = Tree('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_Speciestree_named.newick',format=3)
        spe_tredis = spe_tre.get_distance(mc2node[mc1],mc2node[mc2])/spe_tre.get_farthest_leaf()[1]*100
        print(ko,tre_dis,spe_tredis,n1 in list(mn.traverse()),n2 in list(mn.traverse()),ko2info[ko])

spe_tre = Tree('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_Speciestree_named.newick',format=3)
import itertools
spe_mc2mc_dist = []
for mc1,mc2 in itertools.combinations(list(mc2node),2):
    spe_tredis = spe_tre.get_distance(mc2node[mc1],mc2node[mc2])
    spe_mc2mc_dist.append(spe_tredis)
max_mc_spedist = max(spe_mc2mc_dist)


import numpy as np
ko2sorted_pair = {}
for ko in target_kos:

    #ko = 'K03474'
    tre = ko2tre[ko]
    gene_mc2mc_dist = []
    for mc1,mc2 in itertools.combinations(list(ko2mc2node[ko]),2):
        tre_dis = tre.get_distance(ko2mc2node[ko][mc1],ko2mc2node[ko][mc2])
        gene_mc2mc_dist.append(tre_dis)
    max_mc_genedist = max(gene_mc2mc_dist)
    mc_list = list(ko2mc2node[ko])
    mc2mc_dist = defaultdict(dict)
    for mc1,mc2 in itertools.combinations(mc_list,2):
        if mc1 not in mc2node or mc2 not in mc2node:
            continue
        tre_dis = tre.get_distance(ko2mc2node[ko][mc1],ko2mc2node[ko][mc2])/max_mc_genedist*100
        spe_tredis = spe_tre.get_distance(mc2node[mc1],mc2node[mc2])/max_mc_spedist*100
        mc2mc_dist[mc1][mc2] = (tre_dis,spe_tredis)
    ko2sorted_pair[ko] = sorted([((np.log(v1/v2)),mc1,mc2) for mc1,mcd in mc2mc_dist.items() for mc2,(v1,v2) in mcd.items()])

for k,v in ko2sorted_pair.items():
    print(k,v[:2])

import itertools
for mc1,mc2 in itertools.combinations(mc2node,r=2):
    n1, n2 = mc2node[mc1],mc2node[mc2]
    tre_dis = tre.get_distance(n1,n2)

    for ko,kotre in ko2tre.items():
        largest_ratio = 0
        for n in kotre.traverse():
            all_leafs = n.get_leaf_names()
            genomes = set([_.split('_')[0] for _ in all_leafs if _ in tre.get_leaf_names()])
            if len(genomes)==0:
                continue
            intersec1 = [_ for _ in genomes if _ in n1.get_leaf_names()]
            intersec2 = [_ for _ in genomes if _ in n1.get_leaf_names()]
            if len(intersec)/len(genomes)>=0.7:
                ratio = len(intersec)/len(coral_gs1)
                if ratio >=largest_ratio:
                    tpath2n[ko] = (n,ratio)
                    largest_ratio = ratio
                    ko2inter[ko] = [_ for _ in all_leafs if _.split('_')[0] in intersec]


##### these genes are not clustered in a specific region
genome2info = defaultdict(list)
for ko,(node,ratio) in sorted(tpath2n.items(),key=lambda x:x[1][1]):
    for l in node.get_leaf_names():
        genome = l.split('_')[0]
        if ratio >=0.9:
            genome2info[genome].append((int(l.split('_')[-1].replace('gene','')),l,ko))





