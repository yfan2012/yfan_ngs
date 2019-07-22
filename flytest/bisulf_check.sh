#!/bin/bash

fqdir=/kyber/Data/NGS/Raw/190718_flytest
datadir=/kyber/Data/NGS/projects/190718_flytest

if [ $1 == prep_genome ] ; then
    bismark_genome_preparation --bowtie2 --genomic_composition --parallel 36 $datadir/ref
fi

if [ $1 == bismark ] ; then
    for i in meth unmeth ;
    do
	mkdir -p $datadir/$i
	bismark -o $datadir/$i --prefix $i --parallel 36 -p 36 $datadir/ref -1 $fqdir/$i*R1*.fastq.gz -2 $fqdir/$i*R2*.fastq.gz
    done
fi


if [ $1 == meth_extract ] ; then
    for i in meth unmeth ;
    do
	bismark_methylation_extractor --parallel 36 -p --include_overlap -o $datadir/$i $datadir/$i/$i*.bam
    done
fi

	
