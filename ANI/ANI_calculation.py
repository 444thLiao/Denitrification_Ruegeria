
# ! fastani
from os.path import realpath
import pandas as pd
odir = ''
indir = ''

gids = open('used_genomes.list').read().strip().split('\n')
t = []
for gid in gids:
    fna = f"{indir}/{gid}.fna"
    t.append(realpath(fna))
with open(f'{odir}/genome.list','w') as f1:
    f1.write('\n'.join(t))
cmd = f"fastANI --ql ./genome.list --rl ./genome.list -o pairwise_ANI.out"

a = pd.read_csv(f'{odir}/pairwise_ANI.out',sep='\t',header=None)
a.loc[:,0] = [_.split('/')[-1].replace('.fna','').replace('.fasta','') for _ in a[0]]
a.loc[:,1] = [_.split('/')[-1].replace('.fna','').replace('.fasta','') for _ in a[1]]
# a.groupby([0,1])
ani_df = a.pivot(index=0,columns=1,values=2)
ani_df.to_csv(f'{odir}/445_genome_pairwiseANI.tab',sep='\t',index=1,index_label='genome')
