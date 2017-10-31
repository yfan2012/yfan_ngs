#!/bin/bash

datadir=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample30X/
ref=/mithril/Data/NGS/Reference/cneo/grubii/grubii
dest=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample30X/align/
cov=30

bowtie2 --threads 8 -x ${ref} -1 ${datadir}down_${cov}_C1776_S1_L001_R1_001.fastq -2 ${datadir}down_${cov}_C1776_S1_L001_R2_001.fastq -S ${dest}cneo_1776_${cov}.sam &> ${dest}1776_${cov}_log.txt &
bowtie2 --threads 8 -x ${ref} -1 ${datadir}down_${cov}_C1939_S2_L001_R1_001.fastq -2 ${datadir}down_${cov}_C1939_S2_L001_R2_001.fastq -S ${dest}cneo_1939_${cov}.sam &> ${dest}1939_${cov}_log.txt & 
bowtie2 --threads 8 -x ${ref} -1 ${datadir}down_${cov}_C2898_S3_L001_R1_001.fastq -2 ${datadir}down_${cov}_C2898_S3_L001_R2_001.fastq -S ${dest}cneo_2898_${cov}.sam &> ${dest}2898_${cov}_log.txt


samtools view -bS ${dest}cneo_1776_${cov}.sam > ${dest}cneo_1776_${cov}.bam &
samtools view -bS ${dest}cneo_1939_${cov}.sam > ${dest}cneo_1939_${cov}.bam &
samtools view -bS ${dest}cneo_2898_${cov}.sam > ${dest}cneo_2898_${cov}.bam

samtools sort -o ${dest}cneo_1776_${cov}.sorted.bam ${dest}cneo_1776_${cov}.bam &
samtools sort -o ${dest}cneo_1939_${cov}.sorted.bam ${dest}cneo_1939_${cov}.bam &
samtools sort -o ${dest}cneo_2898_${cov}.sorted.bam ${dest}cneo_2898_${cov}.bam
