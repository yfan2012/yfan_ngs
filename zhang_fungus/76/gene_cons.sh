#!/bin/bash

datadir=~/work/ngs/170928_fungus76
srcdir=~/Code/yfan_ngs/zhang_fungus/76



if [ $1 == rename_fq ] ; then
    python $srcdir/rename.py -i $datadir/fastqs -e R1.fastq.gz -k $datadir/sample_key.csv
    python $srcdir/rename.py -i $datadir/fastqs -e R2.fastq.gz -k $datadir/sample_key.csv
fi

if [ $1 == fix_rename ] ; then
    for i in $datadir/fastqs/*R1.fastq.gz ;
    do
	prefix=`basename $i R1.fastq.gz`
	mv $i $datadir/fastqs/${prefix}_R1.fastq.gz
    done
    for i in $datadir/fastqs/*R2.fastq.gz ;
    do
	prefix=`basename $i R2.fastq.gz`
	mv $i $datadir/fastqs/${prefix}_R2.fastq.gz
    done
fi

if [ $1 == index ] ; then
    for i in $datadir/Reference/ref_alleles/*fasta ;
    do
	prefix=`basename $i .fasta`
	bowtie2-build $i $datadir/Reference/ref_alleles/$prefix
    done
fi

if [ $1 == consensus ] ; then
    touch $datadir/batch_logs
    for i in $datadir/fastqs/*R1.fastq.gz ;
    do
	prefix=`basename $i _R1.fastq.gz`
	sbatch --output=$datadir/batch_logs/$prefix.out --job-name=$prefix $srcdir/consensus.scr $i
    done
fi

if [ $1 == align ] ; then
    mkdir -p $datadir/batch_logs
    for i in $datadir/fastqs/*R1.fastq.gz ;
    do
	prefix=`basename $i _R1.fastq.gz`
	sbatch --output=$datadir/batch_logs/$prefix.align.out --job-name=$prefix $srcdir/align_consensus_parallel.scr $i
    done
fi


if [ $1 == bigmem_realign ] ; then
    for i in $datadir/fastqs/*R1.fastq.gz ;
    do
	prefix=`basename $i _R1.fastq.gz`
	bash $srcdir/align_consensus.scr $i &> $datadir/batch_logs/$prefix.bigmem.out
    done
fi


