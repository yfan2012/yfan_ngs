#!/bin/bash

datadir=/dilithium/Data/NGS/Aligned/170928_fungus76
outdir=~/Dropbox/Lab/fungus_zhang/fungus_76

##This isn't working
parsnp -r /mithril/Data/NGS/Reference/cneo/grubii/grubii.fa -d $datadir/samps_all -p 12 -o $outdir/parsnp_all


##Forgot that you need to delete the newline character at the end of the reference
##for i in grubii jec21 b-3501a ;
for i in jec21 b-3501a ;
do
    ref=/mithril/Data/NGS/Reference/cneo/$i/$i.fa
    parsnp -r $ref -d $datadir/samps_$i -p 12 -o $outdir/parsnp_$i
done

