#!/bin/bash

rawdir=/mithril/Data/NGS/Raw/210510_dunlop_nathan
datadir=/mithril/Data/NGS/projects/dunlop_insert/run2

ref=$datadir/reference/construct.fa

if [ $1 == trim ] ; then
    mkdir -p $datadir/trimmed

    cp ~/software/Trimmomatic-0.39/adapters/all_adapters.fa ~/Code/yfan_ngs/dunlop_insert/run2
    
    for i in NT284 NT285 NT286 NT287;
    do
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 -phred33 \
	     $rawdir/${i}*_L001_R1_001.fastq.gz $rawdir/${i}*_L001_R2_001.fastq.gz \
	     $datadir/trimmed/${i}_fwd_paired.fq.gz $datadir/trimmed/${i}_fwd_unpaired.fq.gz \
	     $datadir/trimmed/${i}_rev_paired.fq.gz $datadir/trimmed/${i}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi

if [ $1 == align ] ; then
    mkdir -p $datadir/align

    bwa index $ref
    
    for i in NT284 NT285 NT286 NT287;
    do
	bwa mem \
	    -t 36 \
	    $ref \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/align/$i.sorted.bam
	samtools index $datadir/align/$i.sorted.bam
    done
fi
