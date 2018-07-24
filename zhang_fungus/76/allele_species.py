import glob
import sys
sys.path.insert(0, '/home/yfan/Code/utils')

from fasta_utils import fasta_dict

'''
take all fastas of alleles as made by cons pipeline, and spit out the allele number that looks most legit
most legit defined as the one with the most non-n characters
'''

consdir='/dilithium/Data/NGS/Aligned/170928_fungus76/cons/cons'
fastas=glob.glob(consdir+'/**/*.fasta')

##read in the species labels
with open('/home/yfan/Dropbox/yfan/fungus_zhang/fungus_76/mlst/species_key.csv') as f:
    species=[x.strip('\r\n').split(',') for x in f.readlines()[1::]]


##reassign species names
for i in species:
    if 'neo' in i[1] or 'grubii' in i[1]:
        i.append('cneo')
    else:
        i.append('gattii')

labels=[x[0] for x in species]
org=[x[2] for x in species]

orglabs=dict(zip(labels, org))

##check for best allele
alleleinfo=[['sample','gene', 'allele', 'len']]
bestalleleinfo=[['sample', 'gene', 'allele', 'len']]
for i in fastas:
    nameseqs=fasta_dict(i)
    gene=i.split('.')[-2]
    samp=i.split('.')[0].split('/')[-1]
    fastaorg=i.split('.')[1].split('_')[0]
    if orglabs[samp] == fastaorg: 
        best=['allele', 0]
        for j in nameseqs:
            alleleinfo.append([samp, gene, j, str(len(nameseqs[j]) - nameseqs[j].count('n'))])
            if len(nameseqs[j]) - nameseqs[j].count('n') > best[1]:
                best[0]=j
                best[1]=len(nameseqs[j]) - nameseqs[j].count('n')
        bestalleleinfo.append([samp, gene, best[0], str(best[1])])



with open('/home/yfan/Dropbox/yfan/fungus_zhang/fungus_76/mlst/mlst_alleles_species.csv', 'w') as f:
    for i in bestalleleinfo:
        if i[-1]!='0':
            f.write(','.join(i)+'\n')



with open('/home/yfan/Dropbox/yfan/fungus_zhang/fungus_76/mlst/mlst_alleles_species_all.csv', 'w') as f:
    for i in alleleinfo:
        if i[-1]!='0':
            f.write(','.join(i)+'\n')

    
            
