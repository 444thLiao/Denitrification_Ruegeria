"""
run at the Hk server

"""

from api_tools.IO_for.itol import get_itoltree
# from evol_tk/bin
from tqdm import tqdm
from collections import defaultdict
from ete3 import Tree
from glob import glob
import pandas as pd
from os.path import *

g2pop = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1783Ruegeria_MCs.tsv',sep='\t',index_col=0)
g2pop = g2pop['MC'].to_dict()
pop2g = defaultdict(list)
for g,pop in g2pop.items():
    pop2g[pop].append(g)



tre = Tree('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/phylogeny/bac120/bac120_MVrooted_2305Ruegeria.newick')

genome_list = pd.read_excel('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_geomes_Feb8.xlsx',index_col=0)
tre.prune(genome_list.index)
tre.write(outfile=f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_Speciestree.newick',format=3)

################# annotated genes
cmds = []
for i in tqdm(glob('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/*.faa')):
    gid = i.split('/')[-1].replace('.faa','')
    cmd = f"/home-user/thliao/software/kofamscan/exec_annotation --cpu 10 -p /mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/annotations/denitrification.hal -k /mnt/home-db/pub/protein_db/kegg/v20230301/ko_list --tmp-dir /mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/annotations/deni_sulfir_out/.{gid}.tmpkofam -o /mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/annotations/deni_sulfir_out/{gid}.kofamout -f mapper-one-line --no-report-unannotated {i}"
    if not exists(f"/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/annotations/deni_sulfir_out/{gid}.kofamout"):
        cmds.append(cmd)

ko2g = pd.read_csv('/mnt/maple/thliao/data/protein_db/kegg/ko_info.tab',
                   sep='\t', header=None, index_col=0)
ko2g.index = [_.split(':')[-1] for _ in ko2g.index]

all_locus2ko = {}
genome2ko = defaultdict(list)
for kofam in tqdm(glob('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/annotations/deni_sulfir_out/*.kofamout')):
    genome = kofam.split('/')[-1].replace('.kofamout','')
    rows = [_.split('\t') for _ in open(kofam).read().strip().split('\n')]
    for locus,ko in rows:
        gene = ko2g.loc[ko,1].split(';')[0].split(',')[0]
        genome2ko[genome].append(gene)
        all_locus2ko[locus]=gene

genome2gene_df = pd.DataFrame(index=tre.get_leaf_names(),columns=['narG','narH','napA','napB','nirS','nirK','norB','norC','nosZ'])

for genome in genome2gene_df.index:
    genome2gene_df.loc[genome,:] = ['Y' if col in genome2ko.get(genome,[]) else 'N' for col in genome2gene_df.columns]


from Bio import SeqIO
genome2gene_df.loc[:,'Genome size/Mbp'] = [round(sum([len(_.seq)/10**6 
for _ in SeqIO.parse(f'/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/ingenomes/{k}.fna','fasta')]),2)
for k in tqdm(genome2gene_df.index)]

genome2gene_df.loc[:,'Ruegeria MC'] = [g2pop[_] for _ in genome2gene_df.index]

allmetadata = pd.read_excel('/mnt/ivy/thliao/project/coral_ruegeria/info_RuegeriaUsed_20240525.xlsx',index_col=0)
fromcoral_idx = allmetadata.index[allmetadata.niche=='coral']
fromcoral_idx = [_ for _ in fromcoral_idx if _ in genome2gene_df.index]
genome2gene_df.loc[fromcoral_idx,'Coral species'] = list(allmetadata.loc[fromcoral_idx,'host name'])
genome2gene_df.loc[fromcoral_idx,'Compartment'] = list(allmetadata.loc[fromcoral_idx,'Compartment'])
genome2gene_df.loc[fromcoral_idx,'Sampling site'] = list(allmetadata.loc[fromcoral_idx,'Location (renamed)'])
genome2gene_df =genome2gene_df.reindex(columns=['Coral species', 'Compartment',
       'Sampling site','Ruegeria MC','Genome size/Mbp',
       'narG', 'narH', 'napA', 'napB', 'nirS', 'nirK', 'norB', 'norC', 'nosZ',
         ])
genome2gene_df.to_excel('../TableS5.xlsx',index=1)

copy_tre = tre.copy()

subdf = genome2gene_df.loc[fromcoral_idx,:]
copy_tre.prune(subdf.index)
copy_tre.write(outfile=f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_Speciestree.newick',format=3)

g2sid = allmetadata['ID'].to_dict()
subdf.loc[:,'strain ID'] = [g2sid[_] for _ in subdf.index]
subdf.to_excel('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/TableS5.xlsx',index=1)


gene2setting = {
'napA': {'shape':'2','color':'#306B34'},
'narG': {'shape':'4','color':'#306B34'},
'nirS': {'shape':'2','color':'#1C5253'},
'nirK': {'shape':'4','color':'#1C5253'},
'norB': {'shape':'2','color':'#F3FFC6'},
'nosZ': {'shape':'2','color':'#C3EB78'},
'complete': {'shape':'2','color':'#CF0013'},}

def justify_complete(genome):
    if set(list(genome2gene_df.loc[genome,['norB','nosZ']])) == set(['Y']):
        if 'Y' in list(genome2gene_df.loc[genome,['napA','narG']]) and 'Y' in list(genome2gene_df.loc[genome,['nirK','nirS']]):
            return True
        else:
            return False       
    return False

from api_tools.itol_func import to_binary_shape
genome2each = defaultdict(list)
for gene,setting in tqdm(gene2setting.items()):
    if gene in genome2gene_df.columns:
        for genome in genome2gene_df.index[genome2gene_df.loc[:,gene]=='Y']:
            genome2each[genome].append(gene)

for g in genome2each:
    if justify_complete(g):
        genome2each[g].append('complete')
text = to_binary_shape(genome2each,info2style=gene2setting,manual_v=['complete',
                                                                     'incomplete'
                                                                     'napA','narG','nirS','nirK','norB','nosZ'],dataset_name=f"denitrification",unfilled_other=True)
with open(f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/denitrification_itolbinary1.txt','w') as f1:
    f1.write(text)
subg2inc = {k:{'incomplete'} for k,v in genome2each.items() if 'complete' not in v}    
text = to_binary_shape(subg2inc,info2style={'incomplete': {'shape':'2','color':'#302A2B'}},manual_v=['incomplete'],no_legend=True,unfilled_other=True)
with open(f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/denitrification_itolbinary2.txt','w') as f1:
    f1.write(text)


cmap_l = {
    "BI": "#E65100",
    "LC": "#535379",
    "SW": "#314f34",
    "YTW": "#811706",
    "NP": "#009bb7",
    "PC" : "#b21900",
    'CI':  '#76FF03',
    'SI': '#a8ca5f',
    'PI': '#b5d2be',
}


# genome2gene_df.loc[genome2gene_df.loc[genome2gene_df['Sampling site']=='Peng Chau',:].index,'Sampling site'] = 'Peng Chau (PC)'
compartment2color = {
'skeleton':'#2f2a2b', 
'tissue':'#566d77', 
'mucus':'#dbdbdb'}
site2color = {
'Yam Tsai Wan (YTW)':"#811706",
'Lo Chau (LC)': "#535379",
'Sham Wan (SW)': "#314f34",
'Ninepin (NP)': "#009bb7",
'Bluff Island (BI)': "#E65100", 
'Peng Chau (PC)': "#b21900"}
coral2color = {'Acropora': '#FFB74D',  
               'Platygyra': '#7b1fa2',
               "Favites": "#03a9f4",
               'Oulastrea': '#00C853',
               }

from api_tools.itol_func import to_color_strip
genome2compartment = {k:v for k,v in subdf['Compartment'].to_dict().items() if v in compartment2color}
text = to_color_strip(genome2compartment,compartment2color,dataset_name='Compartment',other_params={'margin': 50, 'STRIP_WIDTH': 50})
with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification//genome2compartment_colorstrip.txt','w') as f1:
    f1.write(text)
genome2site = {k:v for k,v in subdf['Sampling site'].to_dict().items() if v in site2color}
text = to_color_strip(genome2site,site2color,dataset_name='Sampling site',other_params={'margin': 50, 'STRIP_WIDTH': 50})
with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification//genome2site_colorstrip.txt','w') as f1:
    f1.write(text)

gid2coral_genus = {k:v.split(' ')[0] for k,v in subdf['Coral species'].to_dict().items() if v.split(' ')[0] in coral2color}
text = to_color_strip(gid2coral_genus,
                      coral2color, dataset_name='coral species',
                      other_params={'margin': 50, 'STRIP_WIDTH': 50})
with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/coral_species_colorstrip.txt', 'w') as f1:
    f1.write(text)
    

# for name,cols in tqdm(fulldf.iteritems()):
#     remaining_mcs = [_ for _ in list(cols[~cols.isna()].index) if _.startswith('MC')]
#     subtre = tre.copy()
#     gids = []
#     for mc in remaining_mcs:
#         gids+= pop2g[mc]
#     subtre.prune(gids)
#     treefile = './tmp.newick'
#     subtre.write(outfile=treefile)
#     get_itoltree(treefile,
#                  outfile=f"/mnt/ivy/thliao/project/coral_ruegeria/tmp/{name}.pdf",
#                  anno_files=['/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/genome2MC_lt3_colorstrip.txt',
#                              f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nap.itol.txt',
#                              f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nar.itol.txt',
#                              f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nir.itol.txt',
#                              f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nor.itol.txt',
#                              f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nos.itol.txt',
#                              '/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/renamed_with_MC.txt'
#                              ],
#                  keys={'display_mode':'1',
#                        'tree_x':'500',
#                        'label_display':'1',
#                        'ignore_branch_length':'0',
#                        'line_width':"5"})
#     #break    