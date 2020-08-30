#!/bin/bash

datadir=/kyber/Data/NGS/projects/190513_hardwick
refdir=$datadir/reference
ref=$refdir/CRNE_H99.fasta

if [ $1 == align ] ; then
    ##use bowtie to align to ref
    mkdir -p  $datadir/align

    bowtie2-build \
	-q \
	--threads 32 \
	$ref \
	$refdir/CRNE_H99

    for i in 178 197 1694 6341 ;
    do
	bowtie2 \
	    -p 36 \
	    -x $refdir/CRNE_H99 \
	    -1 $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    -2 $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b |
	    samtools sort -@ 36 -o $datadir/align/$i.sorted.bam

	samtools index $datadir/align/$i.sorted.bam
    done
fi
	    
	    
if [ $1 == freebayes ] ; then
    mkdir -p $datadir/vars

    cd ~/software/freebayes/scripts
    samtools faidx $ref
    
    for i in 178 197 1694 6341 ;
    do
	./freebayes-parallel \
	    <(./fasta_generate_regions.py $ref.fai 100000) 36 \
	    -p 1 \
	    -f $ref \
	    -b $datadir/align/$i.sorted.bam \
	    > $datadir/vars/$i.vcf
    done
fi


if [ $1 == annotsnps ] ; then
    for i in 178 197 1694 6341 ;
    do
	bedtools intersect -wa -wb\
		 -a $datadir/vars/$i.vcf \
		 -b $datadir/reference/CRNE_H99.gff \
		 > $datadir/vars/$i.bed
    done
fi

if [ $1 == snpreport ] ; then
    for i in 178 197 1694 6341 ;
    do
	python ~/Code/utils/snps_report.py \
	       -v $datadir/vars/$i.bed \
	       -o $datadir/vars/$i.csv
    done
fi

if [ $1 == buildsnpeff ] ; then
    cd ~/software/snpEff
    java -jar snpEff.jar build -gff3 -v CRNE_H99
fi
	
if [ $1 == snpeff ] ; then
    for i in 178 197 1694 6341 ;
    do
	java -jar ~/software/snpEff/snpEff.jar \
	     CRNE_H99 \
	     -c ~/software/snpEff/snpEff.config \
	     -download \
	     $datadir/vars/$i.vcf \
	     > $datadir/vars/$i.snpeff.vcf
    done
fi
	 
if [ $1 == snpeff_report ] ; then
    for i in 178 197 1694 6341 ;
    do
	python ~/Code/utils/snpeff_report.py \
	       -v $datadir/vars/$i.snpeff.vcf \
	       -o $datadir/vars/$i.snpeff.csv
    done
fi
	
