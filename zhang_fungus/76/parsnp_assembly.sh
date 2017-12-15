#!/bin/bash

assembledir=~/work/ngs/170928_fungus76/assemble
parsnpdir=~/work/ngs/170928_fungus76/parsnp

mkdir -p $parsnpdir

for i in $assembledir/* ;
do
    prefix=${i#$assembledir/}
    cp $i/contigs.fasta $parsnpdir/$prefix.fasta
done



