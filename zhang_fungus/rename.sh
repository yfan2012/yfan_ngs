#!/bin/bash

datadir=/dilithium/Data/NGS/Aligned/170928_fungus76

for i in $datadir/assemble/* ;
do
    prefix=`echo $i | rev | cut -d / -f 1 | rev`
    ##mv $i/contigs.fasta $i/$prefix.fasta
    cp $i/$prefix.fasta $datadir/parsnp_trimmed/
done


