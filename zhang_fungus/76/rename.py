import argparse
import os

'''take in a dir with *.fasta and change all the names to recognizable sample names'''

parser=argparse.ArgumentParser(description='change apl code to sample name')
parser.add_argument('--indir', '-i', type=str, required=True, help='full path to dir with fastas in it')
parser.add_argument('--ext', '-e', type=str, required=True, help='extension to replace, including the '.' (like .fasta)')
args=parser.parse_args()



##read in the sample key into a dictionary
keyfile='/home/yfan/Dropbox/yfan/fungus_zhang/fungus_76/sample_key.csv'
with open(keyfile, 'r') as f:
    content=f.readlines()

namekey={}
for i in content:
    ##this is just some string manipulation.
    namekey[i.split(',')[0]]=i.split(',')[1].replace(' ', '').split('(')[0]


for i in os.listdir(args.indir):
    if i.endswith(args.ext):
        newname=namekey[i.split('.')[0].split('_')[0]]
        os.rename(args.indir+'/'+i, args.indir+'/'+newname+args.ext)

        
