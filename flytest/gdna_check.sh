#!/bin/bash

datadir=/kyber/Data/NGS/projects/190725_flywtf
fqdir=/kyber/Data/NGS/Raw/190725_flywtf

for i in $datadir/ref/*fna.gz ;
do
    prefix=`basename $i .fna.gz`
    if [ $1 == index ] ; then
	bowtie2-build $i $datadir/ref/$prefix &
	bwa index $i &
    fi

    if [ $1 == align ] ; then
	##align with bowtie
	bwa mem -t 54 $datadir/ref/$prefix.fna.gz $datadir/trimmomatic/meth_fwd_paired.fq.gz $datadir/trimmomatic/meth_rev_paired.fq.gz | samtools view -bS | samtools sort -o $datadir/align/meth.bwa.$prefix.sorted.bam
	bwa mem -t 54 $datadir/ref/$prefix.fna.gz $datadir/trimmomatic/unmeth_fwd_paired.fq.gz $datadir/trimmomatic/unmeth_rev_paired.fq.gz | samtools view -bS | samtools sort -o $datadir/align/unmeth.bwa.$prefix.sorted.bam
	bowtie2 -p 54 -x $datadir/ref/$prefix -1 $datadir/trimmomatic/meth_fwd_paired.fq.gz -2 $datadir/trimmomatic/meth_rev_paired.fq.gz | samtools view -bS | samtools sort -o $datadir/align/meth.bow.$prefix.sorted.bam
	bowtie2 -p 54 -x $datadir/ref/$prefix -1 $datadir/trimmomatic/unmeth_fwd_paired.fq.gz -2 $datadir/trimmomatic/unmeth_rev_paired.fq.gz | samtools view -bS | samtools sort -o $datadir/align/unmeth.bow.$prefix.sorted.bam
    fi
done


if [ $1 == trimmomatic ] ; then
    ##try some v aggressive trimming
    for i in meth unmeth ;
    do
	mkdir -p $datadir/trimmomatic
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 -phred33 \
	     $fqdir/${i}*_R1_*.fastq.gz $fqdir/${i}*_R2_*.fastq.gz \
	     $datadir/trimmomatic/${i}_fwd_paired.fq.gz $datadir/trimmomatic/${i}_fwd_unpaired.fq.gz \
	     $datadir/trimmomatic/${i}_rev_paired.fq.gz $datadir/trimmomatic/${i}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:all.fa:4:10:10 LEADING:32 TRAILING:32 SLIDINGWINDOW:4:30 MINLEN:36 AVGQUAL:32
    done
fi


if [ $1 == fastqc ]; then
    mkdir -p $datadir/fastqc
    fastqc -t 36 -o $datadir/fastqc $datadir/trimmomatic/*fq.gz
fi
	     
