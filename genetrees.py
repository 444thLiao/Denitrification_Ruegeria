
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

g2info = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1783Ruegeria_MCs.tsv',sep='\t',index_col=0)
g2pop = g2info['MC'].to_dict()
pop2g = defaultdict(list)
for g,pop in g2pop.items():
    pop2g[pop].append(g)

genome_list = pd.read_excel('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_geomes_Feb8.xlsx',index_col=0)


ko2g = pd.read_csv('/mnt/maple/thliao/data/protein_db/kegg/ko_info.tab',
                   sep='\t', header=None, index_col=0)
ko2g.index = [_.split(':')[-1] for _ in ko2g.index]

all_locus2ko = {}
genome2ko = defaultdict(list)
for kofam in tqdm(glob('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/annotations/deni_sulfir_out/*.kofamout')):
    genome = kofam.split('/')[-1].replace('.kofamout','')
    if genome not in list(genome_list.index):
        continue
    rows = [_.split('\t') for _ in open(kofam).read().strip().split('\n')]
    for locus,ko in rows:
        gene = ko2g.loc[ko,1].split(';')[0].split(',')[0]
        genome2ko[genome].append((gene,locus))
        all_locus2ko[locus]=gene



from Bio import SeqIO
from glob import glob

### auxilary genes
import os
auxilary_genes = {'napA':['napC','napB','napH','napG','napA','napD','napF'],
                  'norB':['norC','norB','norQ','norD','norE'],
                  'nosZ':['nosL','nosY','nosF','nosD','nosZ','nosR'],
                  'nirS':['nirF','nirC','cobA','COG3794','COG3005','nirS','K07234'],
                  }
for i in ['narG','napA','nirS','nirK','norB','nosZ']:
    os.system(f"mkdir -p ./fullgene_trees/{i}/ref")
    os.system(f"mkdir -p ./fullgene_trees/{i}/tmp")
    os.system(f"mkdir -p ./fullgene_trees/{i}/search")
    os.system(f"mkdir -p ./fullgene_trees/{i}/extract")

ref_locus2gene = {'AB11_00896':'nosL',
                'AB11_00897':'nosY',
                'AB11_00898':'nosF',
                'AB11_00899':'nosD',
                'AB11_00900':'nosZ',
                'AB11_00901':'nosR',
                
                'AB11_00926':'napC',
                'AB11_00927':'napB',
                'AB11_00928':'napH',
                'AB11_00929':'napG',
                'AB11_00930':'napA',
                'AB11_00931':'napD',
                'AB11_00932':'napF',
                
                'AB11_00917':'norC',
                'AB11_00918':'norB',
                'AB11_00919':'norQ',
                'AB11_00920':'norD',
                'AB11_00921':'norE',

                'AB11_00909':'nirF',
                'AB11_00910':'nirC',
                'AB11_00911':'cobA',
                'AB11_00912':'COG3794',
                'AB11_00913':'COG3005',
                'AB11_00914':'nirS',
                'AB11_00915':'K07234',
                }

ref_faa = '/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/AB11.faa'
ref_locus2prot = {ref_locus2gene[_.id]:_ 
                  for _ in SeqIO.parse(ref_faa,'fasta') if _.id in ref_locus2gene}

for keygene,aux in auxilary_genes.items():
    collected_genes = []
    for auxgene in aux:
        collected_genes.append(ref_locus2prot[auxgene])
    SeqIO.write(collected_genes,open(f'./fullgene_trees/{keygene}/ref/{keygene}.faa','w'),'fasta')

# extract target genes from gbk files

flanking_len = 10000
gene2genome2subgbk = defaultdict(dict)
for genome in tqdm(genome_list.index):
    fna = f'/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/ingenomes/{genome}.fna'
    gbk = realpath(fna).replace('.fna','.gbk').replace('.fasta','.gbk')
    if any([genome in v for v in gene2genome2subgbk.values()]):
        continue
    for record in SeqIO.parse(gbk,'genbank'):
        for feature in record.features:
            if feature.type == 'CDS':
                locus = feature.qualifiers['locus_tag'][0]
                if locus in all_locus2ko:
                    gene = all_locus2ko[locus]
                    if gene in ['narG','napA','nirS','nirK','norB','nosZ']:
                        start,end = feature.location.start, feature.location.end
                        extracted_region = record[max(0,start-flanking_len):min(len(record),end+flanking_len)]
                        gene2genome2subgbk[gene][genome] = extracted_region


for gene in gene2genome2subgbk:
    for genome in gene2genome2subgbk[gene]:
        cds_list = []
        for feature in gene2genome2subgbk[gene][genome].features:
            if feature.type == 'CDS':
                cds_list.append((feature.qualifiers['locus_tag'][0],feature.qualifiers['translation'][0]))
        with open(f'./fullgene_trees/{gene}/tmp/{genome}.faa','w') as f1:
            for cds in cds_list:
                f1.write(f">{cds[0]}\n{cds[1]}\n")

# for i in gene2genome2subgbk['nosZ']['O12'].features:
#     if i.type == 'CDS':
#         print(i.qualifiers['locus_tag'][0],i.qualifiers.get('gene',[]))

# perform blast against ref
cmds = []
for gene in ['napA','nirS','norB','nosZ']:
    for faa in glob(f'./fullgene_trees/{gene}/tmp/*.faa'):
        genome = basename(faa).replace('.faa','')
        cmd = f"blastp -query {faa} -subject ./fullgene_trees/{gene}/ref/{gene}.faa -outfmt 6 -out ./fullgene_trees/{gene}/search/{genome}.blast"
        if not exists(f"./fullgene_trees/{gene}/search/{genome}.blast") or os.path.getsize(f"./fullgene_trees/{gene}/search/{genome}.blast") == 0:
            cmds.append(cmd)

from bin.multiple_sbatch import sbatch_all
sbatch_all(cmds,thread_per_tasks=2,prefix_name='blast',batch_size=20)

import pandas as pd
from collections import Counter

subgene2seq = defaultdict(list)
for blastout in glob('./fullgene_trees/*/search/*.blast'):
    gene = blastout.split('/')[-3]
    genome = basename(blastout).replace('.blast','')
    # if gene in ['nirS']:
    #     continue
    with open(blastout) as f1:
        rows = f1.read().strip().split('\n')
    # if len(rows) == 1 and rows[0]=='':
    #     print(blastout)
    #     continue
    blastdf = pd.DataFrame([_.split('\t') for _ in rows],
                           columns=['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
    blastdf.loc[:,'evalue'] = blastdf['evalue'].astype(float)
    blastdf.loc[:,'bitscore'] = blastdf['bitscore'].astype(float)
    blastdf = blastdf.loc[blastdf['evalue'].astype(float)<=1e-5]
    blastdf.sort_values('bitscore',ascending=False,inplace=True)
    blastdf.loc[:,'gene'] = [ref_locus2gene[_] for _ in blastdf['sseqid']]
    headblastdf = blastdf.groupby('qseqid').first()

    
    headblastdf = headblastdf.loc[headblastdf['pident'].astype(float)>=50]
    genecounter = Counter(headblastdf['gene'])
    errors = [k for k,v in genecounter.items() if v!=1]
    headblastdf = headblastdf.loc[~headblastdf['gene'].isin(errors)]

    base_faa = {_.id: _ for _ in SeqIO.parse(f'./fullgene_trees/{gene}/tmp/{genome}.faa','fasta')}
    for _,row in headblastdf.iterrows():
        gene = row['gene']
        locus = row.name
        seq = base_faa[locus]
        subgene2seq[gene].append(seq)


## generate genes within each genome
from collections import defaultdict

genome2keygene2genes = defaultdict(lambda : defaultdict(list))
for gene, seq_list in subgene2seq.items():
    keygene = [k for k,v in auxilary_genes.items() if gene in v][0]
    for seq in seq_list:
        genome = seq.id.split('_')[0]
        genome2keygene2genes[genome][keygene].append(gene)

gene2removed_genomes = defaultdict(list)
for genome, keygene2genes in genome2keygene2genes.items():
    for keygene, genes in keygene2genes.items():
        if len(genes) < len(auxilary_genes[keygene])*0.8:
            gene2removed_genomes[keygene].append(genome)
for keygene, removed_genomes in gene2removed_genomes.items():
    print(keygene,len(removed_genomes))

for gene in subgene2seq:
    keygene = [k for k,v in auxilary_genes.items() if gene in v][0]
    seqs = [_ for _ in subgene2seq[gene] if _.id.split('_')[0] not in gene2removed_genomes[keygene]]
    SeqIO.write(seqs,
                open(f'./fullgene_trees/{keygene}/extract/{gene}.faa','w'),'fasta')



cmds = []
for faa in glob('./fullgene_trees/*/extract/*.faa'):
    gene = basename(faa).replace('.faa','')
    keygene = faa.split('/')[-3]
    cmd = f"mafft --auto {faa} > ./fullgene_trees/{keygene}/extract/{gene}.aln; trimal -in ./fullgene_trees/{keygene}/extract/{gene}.aln -out ./fullgene_trees/{keygene}/extract/{gene}.trim -automated1"
    if not exists(f"./fullgene_trees/{keygene}/extract/{gene}.trim"):
        cmds.append(cmd)

sbatch_all(cmds,thread_per_tasks=2,prefix_name='mafft',batch_size=2)

# concatenate them
with open('./fullgene_trees/genomes.list','w') as f1:
    f1.write('\n'.join(genome_list.index))

for keygene in auxilary_genes:
    cmd = f"python3 /home-user/thliao/script/evol_tk/dating_workflow/bin/concat_aln.py -i ./fullgene_trees/{keygene}/extract/ -o ./fullgene_trees/{keygene}.aln -s trim -ct both -simple -gl ./fullgene_trees/genomes.list"
    os.system(cmd)
    

# B2 ./fullgene_trees/nirS/extract/cobA.trim   cobA is breaked into two CDS.
# BJ3 ./fullgene_trees/nirS/extract/K07234.trim
# B4 ./fullgene_trees/nirS/extract/K07234.trim

os.system(f"rm ./fullgene_trees/*/*.model.gz")
cmds = []
for concataln in glob('./fullgene_trees/*.aln'):
    gene = basename(concataln).replace('.aln','')
    concataln = realpath(concataln)
    cmd = f"iqtree -nt AUTO -m MFP -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -s {concataln} -spp {concataln.replace('.aln','.partition')} -pre {dirname(realpath(concataln))}/{gene}/{gene}"
    cmds.append(cmd)


sbatch_all(cmds,thread_per_tasks=10,prefix_name='iqtree')




# generating itol annotation text files
from api_tools.itol_func import to_binary_shape
from ete3 import Tree
for keygene,aux_genes in auxilary_genes.items():
    tre = Tree(f'./fullgene_trees/{keygene}/{keygene}.contree')

    tip2genes = {}
    for genome,keygene2genes in genome2keygene2genes.items():
        genes = keygene2genes[keygene]
        tip2genes[genome] = genes
    text = to_binary_shape(tip2genes,same_color='#2986cc',dataset_name=f'{keygene} operon',
                           manual_v=aux_genes,unfilled_other=True)
    with open(f'./fullgene_trees/{keygene}/{keygene}.itol_binary.txt','w') as f1:
        f1.write(text)


specialist_genomes = [g for g,pop in g2pop.items() if pop in ['MC1','MC6','MC13','MC87']]
text = to_binary_shape({k:['specialist'] for k in specialist_genomes},
                       info2style={'specialist':{'shape':'3','color':'#ff0000'}},
                       dataset_name='specialist',unfilled_other=True)
with open(f'./fullgene_trees/specialist.itol_binary.txt','w') as f1:
    f1.write(text)

for tre in glob('./fullgene_trees/*/*.contree'):
    os.system(f"FastRoot.py -i {tre} -m MV -o {tre.replace('.contree','.MVrooted.newick')}")