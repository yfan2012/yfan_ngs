'''
Filter fungus allele results based on knowledge about the species
'''

with open('/home/yfan/Dropbox/yfan/fungus_zhang/fungus_76/mlst/species_key_short.csv', 'r') as f:
    species={}
    for line in f:
        (key, val)=line.rstrip("\r\n").split(',')
        species[key] = val

with open('/home/yfan/Dropbox/yfan/fungus_zhang/fungus_76/mlst/mlst_alleles.csv', 'r') as f:
    alleles=[]
    for line in f:
        alleles.append(line.rstrip('\r\n').split(','))


species_alleles=[]
for i in alleles[1::]:
    samp=i[0]
    sampspecies=species[samp]
    if sampspecies in i[1]:
        species_alleles.append(i)
        
with open('/home/yfan/Dropbox/yfan/fungus_zhang/fungus_76/mlst/mlst_alleles_species.csv', 'w') as f:
    for i in species_alleles:
        f.write(','.join(i)+'\n')
