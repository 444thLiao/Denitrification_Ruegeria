# pylint: disable=consider-using-set-comprehension



from ete3 import Tree
import pandas as pd
from collections import defaultdict
from glob import glob
tre = Tree(f'/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/phylogeny/bac120/bac120_MVrooted_2305Ruegeria.newick')
gids = list(tre.get_leaf_names())
gids.remove('FO2')
genome_list = pd.read_excel('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_geomes_Feb8.xlsx',index_col=0)

tre.prune(list(genome_list.index))
gids = tre.get_leaf_names()

g2pop = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1837Ruegeria_MCs.tsv',sep='\t',index_col=0)
g2pop = g2pop['MC'].to_dict()
g2pop = {k:v for k,v in g2pop.items() if k in gids}
pop2g = defaultdict(list)
for g,pop in g2pop.items():
    if g in gids:
        pop2g[pop].append(g)



target_genomes = ['AL3', 'AB11', 'L12', 'AG8', 'A5', 'AN11']


nanopore_sequenceing = [_.split('/')[-1].replace('.faa','') for _ in glob('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/all_faa/*.faa')]

pop2node = defaultdict(list)
nodes = [tre]
while 1:
    if len(nodes) == 0:
        break
    current_node = nodes.pop(0)
    leaves = current_node.get_leaf_names()
    pop = set([g2pop[_] for _ in leaves])
    mc = list(pop)[0]
    if len(pop) ==1:
        pop2node[mc].append(current_node)
    else:
        nodes.extend(current_node.children)

import random
remaining_gids = []
for pop,node_list in pop2node.items():
    for node in node_list:
        leaves = node.get_leaf_names()
        if len(leaves)<=1:
            remaining_gids.extend(leaves)
            continue    
        n1 = [_ for _ in leaves if _ in nanopore_sequenceing]
        if set(leaves).intersection(set(target_genomes)):
            remaining_gids.append(list(set(leaves).intersection(set(target_genomes)))[0])
            continue
        if len(n1) >=1:
            n1 = n1[0]
        else:
            random.shuffle(leaves)
            n1 = leaves[0]
        remaining_gids.append(n1)

print(len(remaining_gids),len(set(remaining_gids)))


prunedtree = tre.copy()
prunedtree.prune(remaining_gids)

import os
from os.path import *
os.system(f'mkdir -p /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/OG_results/infaa/')


faa_list = []
for gid in remaining_gids:
    if gid in nanopore_sequenceing:
        faa_list.append(realpath(f"/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/all_faa/{gid}.faa"))
    else:
        faa_list.append(realpath(f"/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins//{gid}.faa"))
for faa in faa_list:
    os.system(f"ln -sf {faa} /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/OG_results/infaa/")

"orthofinder -og -f ./infaa -o ./of_out -t 30 -a 30"

import time
from glob import glob
while 1:
    og_table = '/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/OG_results/of_out/Results_Feb12/Orthogroups/Orthogroups.tsv'
    if exists(og_table):
        break
    time.sleep(60)
    print('waiting')
og_df = pd.read_csv(og_table,sep='\t',index_col=0)

from Bio import SeqIO
cmds = []
for i in glob(f"{dirname(dirname(og_table))}/Orthogroup_Sequences/*.fa"):
    if len(list(SeqIO.parse(i,'fasta')))<len(remaining_gids)*0.5:
        continue
    cmd = f"mafft {i} > {i.replace('.fa','.aln')}; iqtree -nt 2 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -s {i.replace('.fa','.aln')} -pre {i.replace('.fa','')} "
    cmds.append(cmd)
    
from bin.jobarray_sbatch import generate_sbatch_job_array
generate_sbatch_job_array('job_array.sbatch',
                          glob(f"{dirname(dirname(og_table))}/Orthogroup_Sequences/*.fa"),
                          "mafft {input} > {basename}.aln; iqtree -nt 2 -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -s {basename}.aln -pre {basename} ")





#### ancestral reconstruction  (angst)
angst_odir = '/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/angst'
angst_st = prunedtree.copy()
angst_st.set_outgroup('GNM000011965')
with open(join(angst_odir,'stree.newick'),'w',encoding='utf-8') as f1:
    f1.write('(' + angst_st.write(format=3).replace("NoName",'')[:-1] +':0.);')
infiles = []
for nwk in glob(f'{dirname(angst_odir)}/generax/*.newick'):
    ogname = nwk.split('/')[-1].replace('.newick','')
    if exists(f"{angst_odir}/out/{ogname}/AnGST.newick"):
        continue
    a = open(nwk).read()
    with open(join(angst_odir,f'{ogname}.newick'),'w',encoding='utf-8') as f1:
        f1.write('(' + a.replace("NoName",'')[:-1] +':0.);')
    informat =f'''species={angst_odir}/stree.newick
gene={angst_odir}/{ogname}.newick
penalties={angst_odir}/penalty.file
output={angst_odir}/out/{ogname}'''
    with open(f"{angst_odir}/{ogname}.input",'w',encoding='utf-8') as f1:
        f1.write(informat)
    
        infiles.append(f"{angst_odir}/{ogname}.input")

generate_sbatch_job_array('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/angst_job_array.sbatch',
                          infiles,
                          "/home-user/thliao/anaconda3/envs/py2/bin/python2.7 /home-user/thliao/software/angst/angst_lib/AnGST.py {input}",
                          log_dir='/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/angst/log/',jobname='angst',percpu=2)
#### ancestral reconstruction （ecceTERA）
infiles = []
for ufboot in glob('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/OG_results/of_out/Results_Feb12/Orthogroup_Sequences/*.ufboot'):
    og = ufboot.split('/')[-1].replace('.ufboot','')
    if not exists(f"/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/ecceTERA/{og}/geneTree"):
        infiles.append(ufboot)

cmd = "/home-user/thliao/software/ecceTERA/bin/ecceTERA_linux64 species.file=/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/generax/stree.newick gene.file={input} dated=0 print.reconciliations=2 amalgamate=true sylvx.reconciliation=true output.dir=/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/ecceTERA/{basename} print.newick=true verbose=true print.graph=true print.info=true"
generate_sbatch_job_array('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/ecceTERA_job_array.sbatch',
                          infiles,
                          cmd,
                          log_dir='/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/ecceTERA/logs/',jobname='eccetera',percpu=2)




###########################################
### parse the results
###
###########################################
import os
for f in glob('ecceTERA/OG*/geneTree'):
    if os.path.getsize(f)==0:
        os.system(f"rm {f}")


from ete3 import Tree
from tqdm import tqdm
spec_f = []
for f in tqdm(glob('ecceTERA/OG*/geneTree')):
    try:
        example1 = Tree(f)
    except:
        continue
    add_n = []
    spec = False
    for n in example1.traverse():
        leaves = [_.split('_')[0] for _ in n.get_leaf_names()]
        mc = [g2pop[_] for _ in leaves]
        diffmc = set(mc).difference(set('MC1,MC6,MC13,MC18'.split(',')))
        if len(diffmc)==0 and len(set(mc)) >1:
            spec = True
            add_n.append(n)
        #     print(n,mc)
    if spec:
        spec_f.append((f,add_n))

for (tre,list_n) in spec_f:
    for n in list_n:
        leaves = [_.split('_')[0] for _ in n.get_leaf_names()]
        mc = [g2pop[_] for _ in leaves]
        diffmc = set(mc).intersection(set('MC1,MC6,MC13,MC18'.split(',')))


#### angst

fs = glob('./angst/out/OG*/AnGST.newick')
angst_og = []
for f in tqdm(fs):
    example1 = Tree(f,format=3)
    add_n = []
    spec = False
    for n in example1.traverse():
        leaves = [_.split('_')[0] for _ in n.get_leaf_names()]
        mc = [g2pop[_] for _ in leaves]
        diffmc = set(mc).difference(set('MC1,MC6,MC13,MC18'.split(',')))
        if len(diffmc)==0 and len(set(mc)) >1:
            spec = True
            add_n.append(n)
        #     print(n,mc)
    if spec:
        angst_og.append((f,add_n))
print(len(angst_og))


### 
OG_set1 = [_[0].split('/')[-2] for _ in spec_f ]
OG_set2 = [_[0].split('/')[-2] for _ in angst_og ]
shared_OG = set(OG_set1).intersection(set(OG_set2))
print(len(shared_OG))
union_OG = set(OG_set1).union(set(OG_set2))
print(len(union_OG))


with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/ecceTERA_OG.list','w') as f1:
    f1.write('\n'.join(sorted(set(OG_set1))))
with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/angst_OG.list','w') as f1:
    f1.write('\n'.join(sorted(set(OG_set2))))

kegg_df = pd.read_csv(
    '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno_Revised.tsv', sep='\t', index_col=0)
subkegg_df = kegg_df.loc[[_ for _ in og_df.columns if _ in kegg_df.index],:]
locus2ko = {locus: ko
            for ko, _d in subkegg_df.to_dict().items()
            for genome, locus_list in _d.items()
            for locus in str(locus_list).split(',')}

subog_df = og_df.loc[shared_OG,:]
unionog_df = og_df.loc[union_OG,:]

collect_target_locus = []
ko_collect = []
for og,row in subog_df.iterrows():
    for genome,locus_ in row.to_dict().items():
        if genome in target_genomes and str(locus_) != 'nan':
            collect_target_locus.extend([_.strip() for _ in locus_.split(',')])
        if str(locus_) == 'nan':
            continue
        for l in locus_.split(','):
            if locus2ko.get(l.strip(),'NA')!='NA':
                ko_collect.append(locus2ko.get(l.strip(),'NA'))


eccetera = []
angst = []
for og,row in unionog_df.iterrows():
    for genome,locus_ in row.to_dict().items():
        if genome in target_genomes and str(locus_) != 'nan':
            collect_target_locus.extend([_.strip() for _ in locus_.split(',')])
        if str(locus_) == 'nan':
            continue
        for l in locus_.split(','):
            if locus2ko.get(l.strip(),'NA')!='NA' and og in OG_set1 and og not in OG_set2:
                eccetera.append(locus2ko.get(l.strip(),'NA'))
            if locus2ko.get(l.strip(),'NA')!='NA' and og in OG_set2 and og not in OG_set1:
                angst.append(locus2ko.get(l.strip(),'NA'))

ko2g_df = pd.read_csv(f"/mnt/maple/thliao/data/protein_db/kegg/ko_info.tab",sep='\t',header=None,index_col=0)
ko2g = {ko.split(':')[-1]:str(v).split(';')[0].strip() for ko,v in ko2g_df[1].to_dict().items()}
ko2info = {ko.split(':')[-1]:i for ko,i in ko2g_df[1].to_dict().items()}
ko2info['K26272'] = 'dgcN; D-glutamate N-acetyltransferase [EC:2.3.1.312]'
ko2info['K23997'] = 'nnr; ADP-dependent NAD(P)H-hydrate dehydratase / NAD(P)H-hydrate epimerase [EC:4.2.1.136 5.1.99.6]'





g2ipr_dfs = {g:f.split('/')[-1].split('.')[0] 
           for f in glob('/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/ipr/*.anno/*.tsv')}
_locus2og = {}
allprot = []
for genome,col in tqdm(unionog_df.iteritems()):
    l2faa = {seq.id:seq for seq in SeqIO.parse(f"/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/{genome}.faa",'fasta')}
    for og,locus_ in col.to_dict().items():
        if str(locus_) == 'nan':
            continue        
        for l in locus_.split(','):
            allprot.append(l2faa[l.strip()])
            _locus2og[l] = og
with open('./union_OG.faa','w') as f1:
    SeqIO.write(allprot,f1,'fasta-2line')

cmd = f"export LD_LIBRARY_PATH='' && /home-user/thliao/software/interproscan-5.63-95.0/interproscan.sh -i union_OG.faa -d ./ipr -cpu 20 -pa -goterms -iprlookup -appl CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,SFLD,SMART,SUPERFAMILY,TIGRFAM"

rows = []
for row in open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/AncestralReconstruction/ipr/union_OG.faa.tsv').read().strip().split('\n'):
    rows.append(row.split('\t')) 
databases = set([rows[3] for rows in rows])
db2locus2info = defaultdict(dict)
for row in rows:
    db2locus2info[row[3]][row[0]] = [row[4],row[5],row[11],row[12]]
    db2locus2info['interpro'][row[0]] = [row[11],row[12]]
from collections import Counter
with open('./OG_functions.tsv','w') as f1:
    f1.write(f"OG number\tKEGG number\tKEGG Description\tinterpro number\tinterpro description\tPANTHER number\tPANTHER description\tCDD number\tCDD description\tFound by Software\n")
    for og in union_OG:
        locus_list = [locus.strip() for locus,og_ in _locus2og.items() if og_ == og]
        kos = [locus2ko.get(_,'') for _ in locus_list]
        kos = Counter([_ for _ in kos if _])
        if len(kos)!=0:

            ko = kos.most_common(1)[0][0]
        else:
            ko = 'NA'
        interpro = [db2locus2info['interpro'].get(_,'') for _ in locus_list]
        interpro = [_ for _ in interpro if _ ]
        if len(interpro)!=0:
            p = Counter([_[0] for _ in interpro]).most_common(1)[0][0]
            interpro_info = [_ for _ in interpro if _[0] == p][0]
        else:
            interpro_info = ['NA','NA']

        panther = [db2locus2info['PANTHER'].get(_,'') for _ in locus_list]
        panther = [_ for _ in panther if _ ]
        if len(panther)!=0:
            p = Counter([_[0] for _ in panther]).most_common(1)[0][0]
            panther_info = [_ for _ in panther if _[0] == p][0]
        else:
            panther_info = ['NA','NA']
        cdd = [db2locus2info['CDD'].get(_,'') for _ in locus_list]
        cdd = [_ for _ in cdd if _ ]
        if len(cdd)!=0:
            p = Counter([_[0] for _ in cdd]).most_common(1)[0][0]
            cdd_info = [_ for _ in cdd if _[0] == p][0]
        else:
            cdd_info = ['NA','NA']
        if og in OG_set1 and og in OG_set2:
            t = 'both'
        elif og in OG_set1 and og not in OG_set2:
            t = 'ecceTERA'
        elif og not in OG_set1 and og  in OG_set2:
            t = 'angst'
        f1.write(f"{og}\t{ko}\t{ko2info.get(ko,'NA')}\t{interpro_info[0]}\t{interpro_info[1].capitalize()}\t{panther_info[0]}\t{panther_info[1].capitalize()}\t{cdd_info[0]}\t{cdd_info[1]}\t{t}\n")

d = pd.read_csv('./OG_functions.tsv',sep='\t',index_col=0)
d = d.applymap(lambda x: 'NA' if x == '-' else x)
d.to_excel('OG_functions.xlsx')



##### collect genes and perform interproscan






