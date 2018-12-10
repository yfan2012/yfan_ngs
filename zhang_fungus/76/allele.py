import glob
import sys
##sys.path.insert(0, '/home/yfan/Code/utils')
sys.path.insert(0, '/home-4/yfan7@jhu.edu/Code/utils')
from fasta_utils import fasta_dict

'''
take all fastas of alleles as made by cons pipeline, and spit out the allele number that looks most legit
most legit defined as the one with the most non-n characters
'''

##consdir='/dilithium/Data/NGS/Aligned/170928_fungus76/cons/cons'
consdir='/scratch/groups/wtimp1/ngs/170928_fungus76/batch_cons'
fastas=glob.glob(consdir+'/**/*.fasta', recursive=True)

alleleinfo=[['sample','gene', 'allele', 'len']]
for i in fastas:
    nameseqs=fasta_dict(i)
    best=['allele', 0]
    for j in nameseqs:
        if len(nameseqs[j]) - nameseqs[j].count('n') > best[1]:
            best[0]=j
            best[1]=len(nameseqs[j]) - nameseqs[j].count('n')
    gene=i.split('.')[-2]
    samp=samp=i.split('.')[0].split('/')[-1]
    alleleinfo.append([samp, gene, best[0], str(best[1])])



##with open('/home/yfan/Dropbox/yfan/fungus_zhang/fungus_76/mlst/mlst_alleles.csv', 'w') as f:
with open('/scratch/groups/wtimp1/ngs/170928_fungus76/mlst_alleles.csv', 'w') as f:
    for i in alleleinfo:
        if i[-1]!='0':
            f.write(','.join(i)+'\n')
    
    
            
