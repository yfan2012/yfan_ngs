#!/bin/bash

rawdir=/dilithium/Data/NGS/Raw/170928_NB501653_0049_AHMJ7CBGX2/fastq-fungal_sequencing
datadir=/dilithium/Data/NGS/Aligned/170928_fungus76
refpre=/mithril/Data/NGS/Reference/cneo/b-3501a/b-3501a
ref=$refpre.fasta

for i in $rawdir/*fastq.gz;
do
    prefix=`echo ${i#$rawdir/} | cut -d '_' -f 1,2`
    
    if [ ! -f $datadir/fastqs/${prefix}_R1.fastq.gz ] ; then
	echo catting R1
	cat $rawdir/${prefix}_L001_R1_001.fastq.gz $rawdir/${prefix}_L002_R1_001.fastq.gz $rawdir/${prefix}_L003_R1_001.fastq.gz $rawdir/${prefix}_L004_R1_001.fastq.gz > $datadir/fastqs/${prefix}_R1.fastq.gz ;
    fi

    if [ ! -f $datadir/fastqs/${prefix}_R2.fastq.gz ] ; then
	echo catting R2
	cat $rawdir/${prefix}_L001_R2_001.fastq.gz $rawdir/${prefix}_L002_R2_001.fastq.gz $rawdir/${prefix}_L003_R2_001.fastq.gz $rawdir/${prefix}_L004_R2_001.fastq.gz > $datadir/fastqs/${prefix}_R2.fastq.gz ;
    fi

    if [ ! -f $datadir/align_b-3501a/$prefix.sorted.bam ] ; then
	bowtie2 -p 12 -x $refpre -1 $datadir/fastqs/${prefix}_R1.fastq.gz -2 $datadir/fastqs/${prefix}_R2.fastq.gz | samtools view -bS | samtools sort -o $datadir/align_b-3501a/$prefix.sorted.bam
	samtools index $datadir/align_b-3501a/$prefix.sorted.bam
    fi

    if [ ! -f $datadir/mpileup_b-3501a/$prefix.cons.fq ] ; then
	samtools mpileup -uf $ref $datadir/align_b-3501a/$prefix.sorted.bam | bcftools call -c - | vcfutils.pl vcf2fq > $datadir/mpileup_b-3501a/$prefix.cons.fq
	samtools mpileup -uf $ref $datadir/align_b-3501a/$prefix.sorted.bam | bcftools call -c - | vcfutils.pl varFilter -D100 > $datadir/mpileup_b-3501a/$prefix.vcf
	seqtk seq -a $datadir/mpileup_b-3501a/$prefix.cons.fq > $datadir/cons_b-3501a/$prefix.fasta
    fi
    
done


