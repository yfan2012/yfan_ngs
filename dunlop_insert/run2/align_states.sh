#!/bin/bash

projdir=/mithril/Data/NGS/projects/dunlop_insert
datadir=$projdir/run2


if [ $1 == align ] ; then
    mkdir -p $datadir/check
    bwa index $projdir/refs/ecoli.fa.gz
    
    for i in NT284 NT285 NT286 NT287 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/ecoli.fa.gz \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_bwa_ecoli.sorted.bam
	samtools index $datadir/check/${i}_bwa_ecoli.sorted.bam
    done
fi


if [ $1 == align_withplas ] ; then
    mkdir -p $datadir/check
    bwa index $projdir/refs/ecoli_withplas.fasta.gz
    
    for i in NT284 NT285 NT286 NT287 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/ecoli_withplas.fasta.gz \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_bwa_ecoliplas.sorted.bam
	samtools index $datadir/check/${i}_bwa_ecoliplas.sorted.bam
    done
fi


if [ $1 == align_plas ] ; then
    mkdir -p $datadir/check
    bwa index $projdir/refs/ptkei-dest-cre-loxp-sfgfp.fasta.gz
    
    for i in NT284 NT285 NT286 NT287 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/ptkei-dest-cre-loxp-sfgfp.fasta.gz \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_bwa_plas.sorted.bam
	samtools index $datadir/check/${i}_bwa_plas.sorted.bam
    done
fi


if [ $1 == construct_plas ] ; then
    bwa index $projdir/refs/construct1_plas.fa.gz
    for i in NT284 NT285 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/construct1_plas.fa.gz \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_constructplas.sorted.bam
	samtools index $datadir/check/${i}_constructplas.sorted.bam
    done

    bwa index $projdir/refs/construct2_plas.fa.gz
    for i in NT286 NT287 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/construct2_plas.fa.gz \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_constructplas.sorted.bam
	samtools index $datadir/check/${i}_constructplas.sorted.bam
    done
fi

if [ $1 == ins_plas ] ; then
    bwa index $projdir/refs/ins1_plas.fa
    for i in NT284 NT285 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/ins1_plas.fa \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_insplas.sorted.bam
	samtools index $datadir/check/${i}_insplas.sorted.bam
    done

    bwa index $projdir/refs/ins2_plas.fa
    for i in NT286 NT287 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/ins2_plas.fa \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_insplas.sorted.bam
	samtools index $datadir/check/${i}_insplas.sorted.bam
    done
fi


if [ $1 == count ] ; then
    touch $datadir/run2_alignstates.csv
    echo samp,unmapped,plas_only,ins_only,cre_only,ins_cre,ins_plas,cre_plas,multiple,plas_recomb,cre_recomb,ins_recomb,other_recomb >> $datadir/run2_alignstates.csv
    for i in NT284 NT285 NT286 NT287 ;
    do
	info=`python ~/Code/yfan_ngs/dunlop_insert/check_insert_runs.py \
		            -b $datadir/check/${i}_insplas.sorted.bam`
	echo $info >> $datadir/run2_alignstates.csv
	
    done
    sed -i -e 's/ /\n/g' $datadir/run2_alignstates.csv
fi


if [ $1 == count_ecoli_chr ] ; then
    touch $datadir/run2_ecoli_chr_counts.csv
    for i in NT284 NT285 NT286 NT287 ;
    do
	bam=$datadir/check/${i}_bwa_ecoliplas.sorted.bam
	numreads=`samtools view $bam | awk '$5==60 && $3=="NC_000913.3" {print $0}' | wc -l`
	echo $i,$numreads >> $datadir/run2_ecoli_chr_counts.csv
    done
fi


	
