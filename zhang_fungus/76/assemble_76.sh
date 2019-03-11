#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/fungus/170928_fungus76
fastqdir=$datadir/fastqs



cp ~/software/Trimmomatic-0.38/adapters/NexteraPE-PE.fa ./


for i in $fastqdir/*R1.fastq.gz;
do
    prefix=`basename $i _R1.fastq.gz`
    echo $prefix

    mkdir -p $datadir/assemble/$prefix
    mkdir -p $datadir/assemble/$prefix/trimmed
    
    java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 $i $fastqdir/${prefix}_R2.fastq.gz $datadir/assemble/$prefix/trimmed/${prefix}_forward_paired.fq.gz $datadir/assemble/$prefix/trimmed/${prefix}_forward_unpaired.fq.gz $datadir/assemble/$prefix/trimmed/${prefix}_reverse_paired.fq.gz $datadir/assemble/$prefix/trimmed/${prefix}_reverse_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    spades.py -1 $datadir/assemble/$prefix/trimmed/${prefix}_forward_paired.fq.gz -2 $datadir/assemble/$prefix/trimmed/${prefix}_reverse_paired.fq.gz -t 48 -m 300 -o $datadir/assemble/$prefix
    
done


