#!/bin/bash

datadir=~/work/ngs/170928_fungus76
fastqdir=$datadir/fastqs


mkdir -p $datadir/assemble

for i in $fastqdir/*R1.fastq.gz;
do
    prefix=`echo ${i#$fastqdir/} | cut -d '_' -f 1,2`
    echo $prefix

    mkdir -p $datadir/assemble/$prefix

    spades.py -1 $fastqdir/${prefix}_R1.fastq.gz -2 $fastqdir/${prefix}_R2.fastq.gz -t 70 -m 300 -o $datadir/assemble/$prefix
    
done


