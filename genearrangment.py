
from Bio import SeqIO
from os.path import *

import os
os.chdir('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/keygenestructures')



ref_locus2gene = {
                'AL3_04010':'nosL',
                'AL3_04011':'nosY',
                'AL3_04012':'nosF',
                'AL3_04013':'nosD',
                'AL3_04014':'nosZ',
                'AL3_04015':'nosR',
                
                'AL3_04002':'napC',
                'AL3_04003':'napB',
                'AL3_04004':'napH',
                'AL3_04005':'napG',
                'AL3_04006':'napA',
                'AL3_04007':'napD',
                'AL3_04008':'napF',
                
                'AL3_04031':'norC',
                'AL3_04032':'norB',
                'AL3_04033':'norQ',
                'AL3_04034':'norD',
                'AL3_04035':'norE',

                'AL3_04023':'nirF',
                'AL3_04024':'nirC',
                'AL3_04025':'cobA',
                'AL3_04026':'COG3794',
                'AL3_04027':'COG3005',
                'AL3_04028':'nirS',
                'AL3_04029':'K07234',
                }

flanking = 10000
# MC1-AL3, MC6-AE12, MC13-L12-has Nanopore), non-specialist (MC0-A5, MC15-BG7, MC26 -AN11
info = [('AL3','AL3_2',762531-flanking,797615+flanking,-1),
        ('AE12','contig00004',247503-flanking,279674+flanking,1),
        ('L12','L12_1',3186236-flanking,3217214+flanking,1),
        ('A5','A5_1',258900-flanking,296240+flanking,1),
        ('AN11','AN11_2',855691-flanking,889612+flanking,-1),]
for genome,contig,start,end,DIRECTION in info:
    faa = f'/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/inproteins/{genome}.faa'
    if genome in ['A5',"L12"]:
        faa = f'/home-user/thliao/project/coral_ruegeria/nanopore_processing/canu_o/{genome}/09_prokka/{genome}.faa'
    _gbk = {_.id:_ for _ in SeqIO.parse(realpath(faa).replace('.faa','.gbk'),'genbank')}
    _region = _gbk[contig][start:end]
    if DIRECTION ==-1:
        _region = _region.reverse_complement(id=_region.id,name=_region.name)
        _region.annotations['molecule_type'] =  'DNA'
    with open(f'./{genome}_denitrification.fna','w') as f1:
        SeqIO.write(_region,f1,'fasta-2line')
    with open(f'./{genome}_denitrification.gbk','w') as f1:
        SeqIO.write(_region,f1,'genbank')
    
    with open(f'./{genome}_denitrification.faa','w') as f1:
        for fea in _region.features:
            if fea.type == 'CDS':
                f1.write(f">{fea.qualifiers['locus_tag'][0]}\n{fea.qualifiers['translation'][0]}\n")

full_locus_name = {}
full_locus_name.update(ref_locus2gene)
from collections import Counter
for genome,contig,start,end,DIRECTION in info:
    if genome !='AL3':
        cmd = f"blastp -query ./AL3_denitrification.faa -subject ./{genome}_denitrification.faa -outfmt 6 > AL3_{genome}_blastp.tsv"
    else:
        continue
    #os.system(cmd)
    df = pd.read_csv(f"AL3_{genome}_blastp.tsv",sep='\t',header=None)
    df = df.loc[df[10]<=1e-5]
    df = df.loc[df[0].isin(ref_locus2gene),:]
    df.loc[:,'order'] = [int(_.split('_')[1]) for _ in df[0]]

    df = df.sort_values(11,ascending=False)
    df = df.groupby(0).head(1)
    df = df.sort_values('order')
    # print(len([_ for _ in ref_locus2gene if _ not in list(df[0])]),
    #       len([k for k,v in Counter(df[1]).items() if v>1]))
    for _,row in df.iterrows():
        full_locus_name[row[1]] = ref_locus2gene[row[0]]




_gbk = {_.id:_ for _ in SeqIO.parse(f'Paradenitri.gbk','genbank')}
genomes = [genome for genome,contig,start,end,DIRECTION in info]
# c = []
# for genome1,genome2 in zip(genomes,genomes[1:]):
#     cmd = f"blastn -task blastn -query ./{genome1}_denitrification.fna -subject ./{genome2}_denitrification.fna -outfmt 6 > {genome1}_{genome2}_blastn.tsv"
#     os.system(cmd)
#     df = pd.read_csv(f"{genome1}_{genome2}_blastn.tsv",sep='\t',header=None)
#     #min_identity = df[2].min()
#     for _,row in tqdm(df.iterrows(),total=df.shape[0]):
#         g1, s1,e1 = row[0].split('_')[0],row[6],row[7]
#         g2, s2,e2 = row[1].split('_')[0],row[8],row[9]
#         if abs(e1-s1)>=1000 and row[2]>=30:
#             c.append(row[2])




import pandas as pd
from collections import defaultdict
g2info = pd.read_csv('/mnt/ivy/thliao/project/coral_ruegeria/Merged_popcogeneT/1783Ruegeria_MCs.tsv',sep='\t',index_col=0)
g2pop = g2info['MC'].to_dict()
pop2g = defaultdict(list)
for g,pop in g2pop.items():
    pop2g[pop].append(g)


from pygenomeviz import  Genbank, GenomeViz
from os.path import *
from tqdm import tqdm

gv = GenomeViz(
    fig_width=20,
    fig_track_height=0.7,
    feature_track_ratio=0.3,
    tick_track_ratio=0.2,
    tick_style="axis",
    tick_labelsize=5,
)

gbk_list = [Genbank(f'./{genome}_denitrification.gbk') for genome,contig,start,end,DIRECTION in info ]

def rename_gid(gid):
    gid = gid.split('_')[0]
    if gid in ['AL3','AN11']:
        return f"{gid} ({g2pop[gid]})"
    else:
        return f"{gid} ({g2pop[gid]})"
    
for gbk in gbk_list:
    track = gv.add_feature_track(rename_gid(gbk.name), gbk.range_size, labelsize=15)
    track.add_genbank_features(gbk, plotstyle="arrow",label_type='gene')
    for i in track.features:
        if i.seq_feature.qualifiers.get('locus_tag',['NA'])[0] in full_locus_name:
            i.label = full_locus_name[i.seq_feature.qualifiers['locus_tag'][0]]
            i.labelcolor='red'
            i.facecolor = 'red'
            continue
        i.label = ''

        # label = i.seq_feature.qualifiers.get('gene',['NA'])[0].split('_')[0]
        # if label != 'NA':
        #     i.labelcolor='red'
        #     i.label = label
        #     continue
        # i.label = ''

c = []
for genome1,genome2 in zip(genomes,genomes[1:]):
    cmd = f"blastn -task blastn -query ./{genome1}_denitrification.fna -subject ./{genome2}_denitrification.fna -outfmt 6 > {genome1}_{genome2}_blastn.tsv"
    os.system(cmd)
    df = pd.read_csv(f"{genome1}_{genome2}_blastn.tsv",sep='\t',header=None)
    #min_identity = df[2].min()
    for _,row in tqdm(df.iterrows(),total=df.shape[0]):
        g1, s1,e1 = genome1,row[6],row[7]
        g2, s2,e2 = genome2,row[8],row[9]
        
        if abs(e1-s1)>=1000 and row[2]>=30:
            c.append(row[2])
            gv.add_link((rename_gid(g1), s1, e1), (rename_gid(g2), s2, e2),
                        v=row[2], vmin=50)
fig = gv.plotfig()
_ = gv.set_colorbar(fig,vmin=70,bar_height=0.4,bar_width=0.02)    
_ = fig.savefig('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/keygenestructures/test.pdf')
_ = gv.savefig_html('/mnt/ivy/thliao/project/coral_ruegeria/ipynb/paper_denitrification/keygenestructures/test.html')



