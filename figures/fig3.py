
"""
both protein and dna tree
both full length and amplified region tree

Not include amplicons
"""

from ete3 import Tree
import pandas as pd
from Bio import SeqIO

genome_list = pd.read_excel('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/used_geomes_Feb8.xlsx',index_col=0)
final_genomes = list(genome_list.index)
gene2ko = { 'ATP5B':'K02133',
             'parC':'K02621',
             'nirS2':"K15864"}
kegg_df = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/annotations/KEGG_anno.tab',sep='\t',index_col=0)
missing_gids = [gid for gid in final_genomes if gid not in kegg_df.index]


for gene,ko in gene2ko.items():
    marker_faa = f"/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/{gene}.faa"
    subdf = kegg_df.loc[final_genomes,ko]
    s2s = {_.id:_ for _ in SeqIO.parse(marker_faa,'fasta')}
    print(len([_ for _ in final_genomes if _ not in s2s]))
    break

import os
from subprocess import check_call
from Bio import SeqIO
from os.path import *
from tqdm import tqdm
from ete3 import Tree
from glob import glob
import pandas as pd
from collections import defaultdict
########################################################################################################


rep_otu_files = glob('/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/results/amplicon_*/20230811/Pdada2_output/rep.fasta')

def get_ko(rep_otu,rgene=False):
    ko = rep_otu.split('/')[-5]
    assert 'amplicon_' in ko
    NAME = ko.split('_')[-1]
    if rgene:
        return NAME
    else:
        ko = gene2ko[NAME]
        return ko    
########################################################################################################

# # amplified region
# 118 276 parC BA2_03840 474 bp
# 206 364
# 656 656
# 235 378 pccA GNM900116455_04551 429 bp
# 235 400
# 656 656
# 163 331 ATP5B CH11_02417 504 bp
# 163 331
# 634 634

import pandas as pd
######! protein sequence (trim)
for rep_otu in tqdm(rep_otu_files):
    gene = get_ko(rep_otu,rgene=True)
    if gene !='nirS':
        continue
    pos_fa = join(dirname(rep_otu),'Positive_rep.fasta')
    marker_faa = f"/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/{gene}.faa"
    if not exists(f"{dirname(rep_otu)}/tmp/"):
        os.makedirs(f"{dirname(rep_otu)}/tmp")
    cmd = f"blastx -query {pos_fa} -subject {marker_faa} -outfmt 6 |head -n100 > {dirname(rep_otu)}/tmp/corrected626G_ampliconSearch.tbl"
    #print(cmd)
    #if not exists(f"{dirname(rep_otu)}/tmp/corrected626G_ampliconSearch.tbl"):
    check_call(cmd,shell=1)
    df = pd.read_csv(f'{dirname(rep_otu)}/tmp/corrected626G_ampliconSearch.tbl',sep='\t',header=None)
    subdf = df.loc[df[10]<1e-3,:]
    ref = subdf.iloc[0,:]
    s,e = ref[8],ref[9]
    print(s,e,gene,ref[1], (e-s)*3, 'bp')
    # 288 362 nirS GCF_013031405.1 222 bp
    # 118 276 parC BA2_03840 474 bp
    # 235 393 pccA GNM3892_00090 474 bp
    # 163 331 ATP5B AJ12_01732 504 bp
    aln_prot = {_.id:_ 
                for _ in SeqIO.parse(f'/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/{gene}_prot.aln','fasta')}
    r = aln_prot[ref[1]]
    pos = [_ for _,bp in enumerate(r.seq) if bp !='-']
    new_s = [bp for _,bp in enumerate(r.seq) if bp !='-']
    pos_d = dict(zip(range(len(new_s)),pos))
    ns,ne = pos_d[s],pos_d[e]
    print(ns,ne)
    # 288 487 parC
    # 423 766 pccA
    # 495 782 ATP5B
    with open(f'/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/{gene}_trimmed.prot.aln','w') as f1:
        r = [_[ns:ne] for _ in aln_prot.values()]
        r = [_ for _ in r if len([k for k in str(_.seq) if k!='-'])  > 0.8 * (e-s)]
        print(len(r),len(aln_prot))
        # 817 821 nirS
        # 547 565 parC
        # 549 568 pccA
        # 551 565 pccA
        SeqIO.write(r,f1,'fasta-2line')

    cmd = f"FastTreeMP -lg -gamma /mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/{gene}_trimmed.prot.aln > /mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/phylogeny/{gene}_prot_trimmed_fasttree.newick"
    # if not exists(f"{dirname(ofaa)}/phylogeny/{gene}_Trimmed_656G.contree"):
    check_call(cmd,shell=1)
    
    
    #### rename tree
    gene_tre = Tree(f'/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/phylogeny/{gene}_prot_trimmed_fasttree.newick')
    gene_tre.resolve_polytomy()
    gene_tre.prune([_ for _ in gene_tre.get_leaf_names() if _.split('_')[0] in final_genomes])
    for l in gene_tre.get_leaves():
        l.name = l.name.split('_')[0]
    c = 1
    for n in gene_tre.traverse():
        if not n.is_leaf():
            n.name = f"INode{c}"
            c+=1
    gene_tre.write(outfile=f'/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/phylogeny/{gene}_Trimmed_656G_renamed.newick',format=3)

# rid2seq = {_.id:_ for _ in SeqIO.parse('/mnt/ivy/thliao/project/coral_ruegeria/data_processing/annotations/656G.ffn','fasta')}
# 372 849 parC O6_02136 477 bp
# 607 1084
# 656 656
# 702 1179 pccA O6_00371 477 bp
# 702 1198
# 656 656
# 487 995 ATP5B O6_00119 508 bp
# 487 995
# 634 634

######! nucleotide sequence (remove 3rd)
cmds = []
for aln in glob(f'/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/*_nucl.aln'):
    if 'r3' in aln.split('/')[-1]:continue
    gene = aln.split('/')[-1].split('_')[0] 
    if gene not in gene2ko: continue
    final_odir = f"{dirname(aln)}/"
    ffn_records = {_.id:_ for _ in SeqIO.parse(dirname(aln)+f'/{gene}.ffn','fasta')}
    aln_records = {_.id:_ for _ in SeqIO.parse(aln,'fasta')}
    ridx = [k for k in aln_records if k in ffn_records][0]
    output_aln = final_odir + basename(aln.replace('_nucl.aln',
                                                   '_r3_nucl.aln'))
    seq = remove_3rd_foramplicons(aln_records[ridx],start=0)
    with open(output_aln,'w') as f1:
        seqs = [remove_3rd_foramplicons(s,start=ns) 
                for rid,s in aln_records.items()]
        SeqIO.write(seqs,f1,'fasta-2line')
    final_tree = join(dirname(aln),'phylogeny',
                      basename(output_aln).replace('.aln','.fasttree.newick'))
    #if exists(final_tree):continue
    cmds.append(f"FastTreeMP -gtr -gamma {output_aln} > {final_tree}")   


######! nucleotide sequence (trim)
for rep_otu in tqdm(rep_otu_files):
    gene = get_ko(rep_otu,rgene=True)
    if gene !='nirS':continue
    pos_fa = join(dirname(rep_otu),'Positive_rep.fasta')
    cds_fna = f"/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/{gene}.ffn"
    cmd = f"blastn -query {cds_fna} -subject {pos_fa} -outfmt 6 |head -n100 > {dirname(rep_otu)}/tmp/corrected626G_ampliconSearch_nucl.tbl"
    #print(cmd)
    #if not exists(f"{dirname(rep_otu)}/tmp/corrected626G_ampliconSearch_nucl.tbl"):
    check_call(cmd,shell=1)
    df = pd.read_csv(f'{dirname(rep_otu)}/tmp/corrected626G_ampliconSearch_nucl.tbl',sep='\t',header=None)
    max_v = df[9].max()
    subdf = df.loc[(df[10]<1e-3) & (df[8]==1) & (df[9]==max_v) ,:]
    ref = subdf.iloc[0,:]
    s,e = ref[6],ref[7]
    print(s,e,gene,ref[0], (e-s), 'bp')
    # 860 1086 nirS L7 226 bp
    # 487 995 ATP5B O6 508 bp
    # 372 849 parC O6 477 bp
    # 702 1179 pccA O6 477 bp
    aln_nucl = {_.id:_ 
                for _ in SeqIO.parse(f'{dirname(cds_fna)}/{gene}_nucl.aln','fasta')}
    r = aln_nucl[ref[0]]
    pos = [_ for _,bp in enumerate(r.seq) if bp !='-']
    new_s = [bp for _,bp in enumerate(r.seq) if bp !='-']
    pos_d = dict(zip(range(len(new_s)),pos))
    ns,ne = pos_d[s],pos_d[e]
    print(ns,ne,gene)
    # 1527 1774 nirS
    # 288 487 parC
    # 423 766 pccA
    # 495 782 ATP5B
    with open(f'{dirname(cds_fna)}/{gene}_trimmed.nucl.aln','w') as f1:
        r = [_[ns:ne] for _ in aln_nucl.values()]
        r = [_ for _ in r if len([k for k in str(_.seq) if k!='-'])  > 0.8 * (e-s)]
        print(len(r),len(aln_nucl))
        # 807 808 nirS
        # 547 565 parC
        # 549 568 pccA
        # 551 565 pccA
        SeqIO.write(r,f1,'fasta-2line')

from ..tk import *
######! nucleotide sequence (trim + remove 3rd)
for aln in glob(f'/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/*_trimmed.nucl.aln'):
    gene = aln.split('/')[-1].split('_')[0] 
    final_odir = f"{dirname(aln)}/"
    faa_records = {_.id:_ for _ in SeqIO.parse(dirname(aln)+f'/{gene}.faa','fasta')}
    ffn_records = {_.id:_ for _ in SeqIO.parse(dirname(aln)+f'/{gene}.ffn','fasta')}
    trimmed = {_.id:_ for _ in SeqIO.parse(aln,'fasta')}
    r = [k for k in trimmed if k in ffn_records][0]
    i = 0
    for bp in trimmed[r].seq:
        if bp =='-':
            i+=1
        else:
            break    
    ns,ne,s,e =  map_two_aln(trimmed[r],ffn_records[r])
    ns = ns-i
    output_aln = final_odir + basename(aln.replace('_trimmed.nucl.aln',
                                                   '_trimmed.r3.nucl.aln'))
    seq = remove_3rd_foramplicons(trimmed[r],start=ns)
    with open(output_aln,'w') as f1:
        seqs = [remove_3rd_foramplicons(s,start=ns) 
                for rid,s in trimmed.items()]
        SeqIO.write(seqs,f1,'fasta-2line')
    final_tree = join(dirname(aln),'phylogeny',
                      basename(output_aln).replace('.aln','.fasttree.newick'))
    #if exists(final_tree):continue
    os.system(f"FastTreeMP -gtr -gamma {output_aln} > {final_tree}")   


##! batch iqtree
cmds1 = []
cmds = []
for aln in glob(f'/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/*prot.aln'):
    name = basename(aln).replace('.aln','')
    if not exists(f"{dirname(aln)}/iqtree/{name}.contree"):
        cmds.append(f"iqtree -nt 10 -m TESTMERGE -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {dirname(aln)}/iqtree/{name} -s {aln} ")
        
    else:
        print(name)
    final_tree = join(dirname(aln),'phylogeny',
                      basename(aln).replace('.aln','.fasttree.newick'))  
    if not exists(final_tree):  
        cmds1.append(f"FastTreeMP -lg -gamma {aln} > {final_tree}")
# cmds = []
for aln in glob(f'/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/*nucl.aln'):
    name = basename(aln).replace('.aln','')
    if not exists(f"{dirname(aln)}/iqtree/{name}.contree"):
        cmds.append(f"iqtree -nt 10 -m MFP -mset GTR,HKY85,K80 -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -pre {dirname(aln)}/iqtree/{name} -s {aln} ")
    else:
        print(name)
    final_tree = join(dirname(aln),'phylogeny',
                      basename(aln).replace('.aln','.fasttree.newick'))  
    if not exists(final_tree):  
        cmds1.append(f"FastTreeMP -gtr -gamma {aln} > {final_tree}")
for _ in cmds:
    os.system(_)
    
sbatch_all(cmds,thread_per_tasks=10,prefix_name='iqtree')



for tpath in sorted(glob('/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/iqtree/*.contree')):
    tre = Tree(tpath)
    mp_n = tre.get_midpoint_outgroup()
    tre.set_outgroup(mp_n)      
    named = renamed_tree(tre)
    named.write(outfile=tpath.replace('.contree','.formatted.newick'),
                format=3)

from api_tools.IO_for.itol import get_itoltree
for tpath in sorted(glob('/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/prepared_trees/*.formatted.newick')):
    get_itoltree(tpath,
                basename(tpath.replace('.formatted.newick','.png')),
                anno_files=['/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/genera_colorstrip.txt',
                            '/mnt/ivy/thliao/project/coral_ruegeria/data_processing/phylogeny/itol/popcogent_paperID_lt3.txt'],
                keys = {'line_width':3})






id2seq = {_.id:_ for _ in SeqIO.parse('/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/ATP5B_trimmed.nucl.aln','fasta')}
print(len([_ for _ in id2seq['B4'].seq if _!='']))
# 508 bp
with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/ATP5B_amplicon.nucl.aln','w') as f1:
    SeqIO.write([id2seq[_] for _ in gids if _ in id2seq],f1,'fasta-2line')

id2seq = {_.id:_ for _ in SeqIO.parse('/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/markers_trees/parC_trimmed.nucl.aln','fasta')}
print(len([_ for _ in id2seq['B4'].seq if _!='']))
# 486 bp
with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/parC_amplicon.nucl.aln','w') as f1:
    SeqIO.write([id2seq[_] for _ in gids if _ in id2seq],f1,'fasta-2line')



# nirS2 

seqs = []
for i in tqdm(final_genomes):
    if not exists(f"/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/annotations/deni_sulfir_out/{i}.kofamout"):
        print(i)
    locus2ko = {row.split('\t')[0]:row.split('\t')[1] 
                for row in open(f"/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/annotations/deni_sulfir_out/{i}.kofamout").read().split('\n') if row}
    nirS_locus = [k for k,v in locus2ko.items() if v =='K15864']
    if len(nirS_locus)==1:
        faa = {_.id:_ for _ in SeqIO.parse(f'/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/{i}.faa','fasta')}
        seqs.append(faa[nirS_locus[0]])
    elif len(nirS_locus)>1:
        print(i,'have two nirS')
faa = {_.id:_ for _ in SeqIO.parse(f'/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/BJ3.faa','fasta')}
seqs.append(faa['BJ3_03836'])
with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nirS_prot.faa','w') as f1:
    SeqIO.write(seqs,f1,'fasta-2line')

nucl_seqs = []
for l in seqs:
    locus = l.id
    gid = locus.split('_')[0]
    ffn = {_.id:_ for _ in SeqIO.parse(realpath(f'/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/{gid}.faa').replace('.faa','.ffn'),'fasta')}
    nucl_seqs.append(ffn[locus])
with open('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nirS_nucl.ffn','w') as f1:
    SeqIO.write(nucl_seqs,f1,'fasta-2line')

#####
marker_faa = f"/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nirS_prot.faa"
rep_otu = '/mnt/ivy/thliao/project/coral_ruegeria/data_otheramplicons/results/amplicon_nirS2/20230811/Pdada2_output/rep.fasta'
pos_fa = join(dirname(rep_otu),'Positive_rep.fasta')

if not exists(f"{dirname(rep_otu)}/tmp/"):
    os.makedirs(f"{dirname(rep_otu)}/tmp")
cmd = f"blastx -query {pos_fa} -subject {marker_faa} -outfmt 6 |head -n100 > {dirname(rep_otu)}/tmp/corrected626G_ampliconSearch.tbl"
check_call(cmd,shell=1)

df = pd.read_csv(f'{dirname(rep_otu)}/tmp/corrected626G_ampliconSearch.tbl',sep='\t',header=None)
subdf = df.loc[df[10]<1e-3,:]
ref = subdf.iloc[0,:]
s,e = ref[8],ref[9]
print(s,e,gene,ref[1], (e-s)*3, 'bp')
# 112 223 nirS2 GCF_017872595.1 333 bp
aln_prot = {_.id:_ 
            for _ in SeqIO.parse(f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nirS_prot.aln','fasta')}
r = aln_prot[ref[1]]
pos = [_ for _,bp in enumerate(r.seq) if bp !='-']
new_s = [bp for _,bp in enumerate(r.seq) if bp !='-']
pos_d = dict(zip(range(len(new_s)),pos))
ns,ne = pos_d[s],pos_d[e]
print(ns,ne)
# 193 304 nirS2
with open(f'/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nirS2_trimmed.prot.aln','w') as f1:
    r = [_[ns:ne] for _ in aln_prot.values()]
    r = [_ for _ in r if len([k for k in str(_.seq) if k!='-'])  > 0.8 * (e-s)]
    print(len(r),len(aln_prot))
    # 381 381 nirS2
    SeqIO.write(r,f1,'fasta-2line')



pos_fa = join(dirname(rep_otu),'Positive_rep.fasta')
cds_fna = f"/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nirS_nucl.ffn"
cmd = f"blastn -query {cds_fna} -subject {pos_fa} -outfmt 6 |head -n100 > {dirname(rep_otu)}/tmp/corrected626G_ampliconSearch_nucl.tbl"
check_call(cmd,shell=1)
df = pd.read_csv(f'{dirname(rep_otu)}/tmp/corrected626G_ampliconSearch_nucl.tbl',sep='\t',header=None)
max_v = df[9].max()
subdf = df.loc[(df[10]<1e-3) & (df[8]==1) & (df[9]==max_v) ,:]
ref = subdf.iloc[0,:]
s,e = ref[6],ref[7]
print(s,e,gene,ref[0], (e-s), 'bp')
# 332 669 nirS2 A5_03172 337 bp
aln_nucl = {_.id:_ 
            for _ in SeqIO.parse(f'{dirname(cds_fna)}/nirS_nucl.aln','fasta')}
r = aln_nucl[ref[0]]
pos = [_ for _,bp in enumerate(r.seq) if bp !='-']
new_s = [bp for _,bp in enumerate(r.seq) if bp !='-']
pos_d = dict(zip(range(len(new_s)),pos))
ns,ne = pos_d[s],pos_d[e]
print(ns,ne,gene)
# 549 886 nirS2
with open(f'{dirname(cds_fna)}/{gene}_trimmed.nucl.aln','w') as f1:
    r = [_[ns:ne] for _ in aln_nucl.values()]
    r = [_ for _ in r if len([k for k in str(_.seq) if k!='-'])  > 0.8 * (e-s)]
    print(len(r),len(aln_nucl))
    # 381 381 nirS
    r = [_ for _ in r if _.id.split('_')[0]!='GNM000011965']
    for _ in r:
        _.id = _.id.split('_')[0]
    SeqIO.write(r,f1,'fasta-2line')

cmd = f"iqtree -nt 15 -m MFP -redo -mset GTR,HKY85,K80 -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -s /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nirS2_trimmed.nucl.aln -pre /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nirS2_amplicon.nucl_iqtree"
cmd = 'FastRoot.py -i /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nirS2_amplicon.nucl_iqtree.contree -m MV -o /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/nirS2_amplicon.nucl.MV.newick'
cmd = f"iqtree -nt 15 -m MFP -redo -mset GTR,HKY85,K80 -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -s /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/parC_amplicon.nucl.aln -pre /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/parC_amplicon.nucl_iqtree"


cmd = f"iqtree -nt 15 -m MFP -redo -mset GTR,HKY85,K80 -mrate E,I,G,I+G -mfreq FU -wbtl -bb 1000 -s /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/ATP5B_amplicon.nucl.aln -pre /mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/ATP5B_amplicon.nucl_iqtree"




