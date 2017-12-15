#!/bin/bash

assembledir=/dilithium/Data/NGS/Aligned/170928_fungus76/assemble
parsnpdir=/dilithium/Data/NGS/Aligned/170928_fungus76/assemble/parsnp_assemblies

mkdir -p $parsnpdir

##move all assemblies into one folder
for i in $assembledir/APL* ;
do
    prefix=${i#$assembledir/}
    cp $i/contigs.fasta $parsnpdir/$prefix.fasta
done


##rename so parsnp output is understandable
python ~/Code/yfan_ngs/zhang_fungus/76/rename.py -i $parsnpdir


ref=/mithril/Data/NGS/Reference/cneo/grubii/grubii.fa
outdir=~/Dropbox/Lab/fungus_zhang/parsnp_assemblies

mkdir -p $outdir

parsnp -r $ref -d $parsnpdir -p 12 -o $outdir -c

