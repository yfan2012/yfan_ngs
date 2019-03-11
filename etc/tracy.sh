#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/tracy

mkdir -p $datadir/assemblies

ml python/2.7

for i in $datadir/*fastq.gz ;
do
    prefix=` basename $i .fastq.gz `

    if [ $1 == assemble ] ; then
	spades.py -s $i -t 36 -m 300 -o $datadir/assemble/$prefix
    fi
    if [ $1 == rename ] ; then
	mv $datadir/assemble/$prefix/contigs.fasta $datadir/assemble/$prefix/$prefix.fasta
    fi

    if [ $1 == fix ] && [ $prefix == 2048pa_S35_L001_R1_001 ] ; then
	spades.py -s $i -t 36 -m 300 -o $datadir/assemble/$prefix
    fi
done

