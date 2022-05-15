#!/bin/bash

rawdir=/mithril/Data/NGS/Raw/220418_dunlop_nathan
projdir=/mithril/Data/NGS/projects/dunlop_insert
datadir=$projdir/run5

mkdir -p $datadir

ref1=$projdir/refs/construct1.fa ##magnets
ref2=$projdir/refs/construct2.fa ##ilid
plasref=$projdir/refs/ptkei-dest-cre-loxp-sfgfp.fasta

if [ $1 == trim ] ; then
    mkdir -p $datadir/trimmed

    cp ~/software/Trimmomatic-0.39/adapters/all_adapters.fa ~/Code/yfan_ngs/dunlop_insert/run4
    
    for i in nt398 nt399 nt400 nt401 nt402 nt403 ;
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

    bwa index $ref1
    

    for i in nt398 nt399 nt400 nt401 nt402 nt403 ;
    do
	bwa mem \
	    -t 36 \
	    $ref1 \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/align/$i.sorted.bam
	samtools index $datadir/align/$i.sorted.bam
    done

fi



