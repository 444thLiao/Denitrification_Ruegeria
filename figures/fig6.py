
import pandas as pd
from collections import defaultdict


genomic_component_file = '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/DFD/s2out/100Genomes_pos.tsv'
all_genomic_component_df = pd.read_csv(genomic_component_file,sep='\t',index_col=0)

kegg_df = pd.read_csv(
    '/mnt/ivy/thliao/project/coral_ruegeria/nanopore_processing/annotations/KEGG_anno_Revised.tsv', sep='\t', index_col=0)
g2pop = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1837Ruegeria_MCs.tsv',sep='\t',index_col=0)
g2pop = g2pop['MC'].to_dict()
pop2g = defaultdict(list)
for g,pop in g2pop.items():
    pop2g[pop].append(g)
g2chrom = defaultdict(list)
for row in open('/home-user/thliao/project/coral_ruegeria/nanopore_processing/canu_o/chromosome.txt').read().strip().split('\n'):
    genome,contig = row.split('\t')
    g2chrom[genome].append(contig)

target_genes = {'K00370': 'narG',
                'K00371': 'narH',
                'K02567': 'napA',
                'K02568': 'napB',
                'K00368': 'nirK',
                'K15864': 'nirS',
                'K04561': 'norB',
                'K02305': 'norC',
                'K00376': 'nosZ',
                'K00372': 'nasA',
                # 'K00362': 'nirB',
                # 'K00360': 'nasB'
                }


subdf = kegg_df.loc[:,target_genes]

subdf.columns = [target_genes.get(_,_) for _ in subdf.columns]

def trans(locus):
    if locus not in all_genomic_component_df.index and ',' not in str(locus):
        return locus
    genome = locus.split('_')[0]
    r = []
    for l in locus.split(','):
        row = all_genomic_component_df.loc[l,:]
        if row['contig'] in g2chrom[genome]:
            r.append(f"{l} (chrom)")
        else:
            r.append(f"{l} (plasmid)")
    return ' '.join(r)
    

subdf = subdf.applymap(trans)
subdf.loc[:,'MC'] = [g2pop.get(_,'') for _ in subdf.index]
subdf = subdf.loc[subdf.MC!='',:].fillna('')
subdf.loc[subdf['MC'].isin(['MC13','MC6'])]
subdf.to_excel('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/Denitrification_genes_dis.xlsx')
#################################### 
from os.path import *
from Bio import SeqIO
genedf = kegg_df.loc[:,target_genes]
genedf.columns = [target_genes.get(_,_) for _ in genedf.columns]
gene2seq = defaultdict(list)
for gid,row in genedf.iterrows():
    ffn = realpath(f"/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/ingenomes/{gid}.fna").replace('.fna','.ffn')
    if not exists(ffn):
        continue
    id2seq = {_.id:_ for _ in SeqIO.parse(ffn,'fasta')}
    for gene,v in row.to_dict().items():
        if v in id2seq:
            seq = id2seq[v]
            seq.id = gid
            gene2seq[gene].append(seq)
import os
for gene,seqs in gene2seq.items():
    with open(f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/gene_test/{gene}.fna','w') as f1:
        SeqIO.write(seqs,f1,'fasta-2line')
    os.system(f"mafft /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/gene_test/{gene}.fna > /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/gene_test/{gene}.aln")



from os.path import *
from Bio import SeqIO
genedf = kegg_df.loc[:,target_genes]
genedf.columns = [target_genes.get(_,_) for _ in genedf.columns]
gene2seq = defaultdict(list)
for gid,row in genedf.iterrows():
    ffn = realpath(f"/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/ingenomes/{gid}.fna").replace('.fna','.faa')
    if not exists(ffn):
        continue
    id2seq = {_.id:_ for _ in SeqIO.parse(ffn,'fasta')}
    for gene,v in row.to_dict().items():
        if v in id2seq:
            seq = id2seq[v]
            seq.id = gid
            gene2seq[gene].append(seq)
import os
for gene,seqs in gene2seq.items():
    with open(f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/gene_test/{gene}.faa','w') as f1:
        SeqIO.write(seqs,f1,'fasta-2line')
    os.system(f"mafft /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/gene_test/{gene}.faa > /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/gene_test/{gene}_prot.aln")



def compare_seq(s1,s2):
    total_num = 0
    diff_num = 0
    for b1,b2 in zip(s1.seq,s2.seq):
        if b1 !='-' and b2 !='-':
            total_num+=1
            if b1 !=b2:
                diff_num +=1
    return diff_num,total_num

gene = 'napA'
aln_seqs = {_.id:_ for _ in SeqIO.parse(f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/gene_test/{gene}.aln','fasta')}
print(compare_seq(aln_seqs['AB11'],aln_seqs['L12']))
aln_seqs = {_.id:_ for _ in SeqIO.parse(f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/gene_test/{gene}_prot.aln','fasta')}
print(compare_seq(aln_seqs['AB11'],aln_seqs['L12']))











