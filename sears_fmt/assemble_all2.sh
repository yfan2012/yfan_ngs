#!/bin/bash

##set up some directory names
datadir=~/work/ngs/171115_sears_fmt
fastqdir=$datadir/fastqs

mkdir -p $datadir/assemble


##do the following commands for all the samples in the directory
for i in $fastqdir/*R1_001.fastq.gz;
do
    ##get the sample name from the file of reads
    prefix=`echo ${i#$fastqdir/} | cut -d '_' -f 1,2`
    echo $prefix

    mkdir -p $datadir/assemble/$prefix

    ##assemble using spades
    spades.py -1 $fastqdir/${prefix}_L001_R1_001.fastq.gz -2 $fastqdir/${prefix}_L001_R2_001.fastq.gz -t 70 -m 300 -o $datadir/assemble/$prefix
    
done
