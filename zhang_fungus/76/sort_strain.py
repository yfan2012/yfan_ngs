import os
from shutil import copyfile

datadir='/dilithium/Data/NGS/Aligned/170928_fungus76/'

statfiles=os.listdir(datadir + 'align_grubii')

cneo=['grubii', 'jec21', 'b-3501a']
samplist=[x.split('.')[0] for x in statfiles]
samplist=set(samplist)

for i in samplist:
    ##alignpercents
    align_perc={}
    for strain in cneo:
        with open(datadir + 'stat_' + strain + '/' + i + '.stat.txt') as f :
            content=f.readlines()
        align_perc[strain]=(float(content[4].split(' ')[0])/float(content[0].split(' ')[0]))

    ##https://stackoverflow.com/questions/36502505/keyword-functions-for-python-min-max
    highest=max(align_perc, key=align_perc.get)
    copyfile(datadir+'mpileup_'+highest+'/'+i+'.fasta', datadir+'samps_'+highest+'/'+i+'.fasta')

    print '\t'.join([i, highest, str(align_perc[highest])])


        
    
