#!/bin/bash

##Analysis of grubii variant.
ref=/mithril/Data/NGS/Reference/cneo/grubii/grubii.fa
aligned=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample30X/align/
destdir=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample30X/mpileup/


if [ 1 -eq 0 ]; then
##samtools faidx done on ref
samtools mpileup -uf $ref ${aligned}cneo_1776_50.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${destdir}cneo_1776_50.cons.fq &
samtools mpileup -uf $ref ${aligned}cneo_1939_50.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${destdir}cneo_1939_50.cons.fq &
samtools mpileup -uf $ref ${aligned}cneo_2898_50.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${destdir}cneo_2898_50.cons.fq &

																	 
samtools mpileup -uf $ref ${aligned}cneo_1776_50.sorted.bam | bcftools call -c - | vcfutils.pl varFilter -D100 > ${destdir}cneo_1776_50.vcf &
samtools mpileup -uf $ref ${aligned}cneo_1939_50.sorted.bam | bcftools call -c - | vcfutils.pl varFilter -D100 > ${destdir}cneo_1939_50.vcf &
samtools mpileup -uf $ref ${aligned}cneo_2898_50.sorted.bam | bcftools call -c - | vcfutils.pl varFilter -D100 > ${destdir}cneo_2898_50.vcf
fi


##convert to fasta to appease parsnp
mkdir ${destdir}fastas
seqtk seq -a ${destdir}cneo_1776_50.cons.fq > ${destdir}fastas/cneo_1776_50.fasta
seqtk seq -a ${destdir}cneo_1939_50.cons.fq > ${destdir}fastas/cneo_1939_50.fasta
seqtk seq -a ${destdir}cneo_2898_50.cons.fq > ${destdir}fastas/cneo_2898_50.fasta
