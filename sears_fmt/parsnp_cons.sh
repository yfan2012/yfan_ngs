#!/bin/bash

datadir=/dilithium/Data/NGS/Aligned/171115_sears_fmt
statdir=$datadir/align_stats

mkdir -p $statdir

if [ $1 == old ] ; then
    ##Check the alignment rates. I think I saw some really low ones.
    for i in $datadir/*.sorted.bam ;
    do
	prefix=`echo ${i%.sorted.bam} | rev | cut -d '/' -f 1 | rev`
	samtools flagstat $i > $statdir/$prefix.stat.txt
    done
    
    
    consdir=$datadir/mpileup
    mkdir -p $consdir
    
    cp $datadir/*fasta $consdir/
    
    
    ref=/mithril/Data/NGS/Reference/ecoli/ecoli.fasta
    outdir=~/Dropbox/Lab/sears_fmt/parsnp_cons
    mkdir -p $outdir
    
    
    parsnp -r $ref -d $consdir -p 12 -o $outdir -c
fi


    
