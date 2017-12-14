#!/bin/bash

rawdir=/dilithium/Data/NGS/Raw/171115_sears_fmt
datadir=/dilithium/Data/NGS/Aligned/171115_sears_fmt
refpre=/mithril/Data/NGS/Reference/ecoli/ecoli


for i in $rawdir/*R1_001.fastq.gz;
do
    ##
    prefix=`echo ${i#$rawdir/} | cut -d '_' -f 1,2`

    ##Read mapping/alignment
    if [ ! -f $datadir/align_grubii/$prefix.sorted.bam ] ; then
	bowtie2 -p 12 -x $refpre -1 $datadir/fastqs/${prefix}_R1.fastq.gz -2 $datadir/fastqs/${prefix}_R2.fastq.gz | samtools view -bS | samtools sort -o $datadir/align_grubii/$prefix.sorted.bam
	samtools index $datadir/align_grubii/$prefix.sorted.bam
    fi

    
    if [ ! -f $datadir/mpileup_grubii/$prefix.cons.fq ] ; then
	samtools mpileup -uf $ref $datadir/align_grubii/$prefix.sorted.bam | bcftools call -c - | vcfutils.pl vcf2fq > $datadir/mpileup_grubii/$prefix.cons.fq
	samtools mpileup -uf $ref $datadir/align_grubii/$prefix.sorted.bam | bcftools call -c - | vcfutils.pl varFilter -D100 > $datadir/mpileup_grubii/$prefix.vcf
	seqtk seq -a $datadir/mpileup_grubii/$prefix.cons.fq > $datadir/cons_grubii/$prefix.fasta
    fi
    
done


