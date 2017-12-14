#!/bin/bash

##Downsample reads from a zipped fastq file
##For 40X coverage of 19M genome on 150bp long paired reads: 25000000ish paired end reads
cov=10

#number of reads you want to have per file at the end (calculate based on coverage)
num=650000

##List of files
data=/atium/Data/NGS/Raw/160921_SSSSS_Crypt/C*.fastq.gz

#Destination directoryp
downdest=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample${cov}X/
mkdir ${downdest}

#File name modifier for downsampled fastq
modif=down_${cov}_


echo Downsampling---------------------------------
for i in ${data}
do
    #Get rid of path stuff
    file=${i##*/}

    #Get rid of gz file extension
    name=${file%.*}

    #seed
    seed=${file:1:3}
    
    #Downsample
    #-s is seed. Must be same for read pairs. I think it can be the same for all. 
    seqtk sample -s100 ${i} ${num} > ${downdest}/${modif}${name}
done

echo Downsampling done----------------------------



echo Aligning-------------------------------------
ref=/mithril/Data/NGS/Reference/cneo/grubii/grubii
aligndest=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample${cov}X/align/
mkdir ${aligndest}

bowtie2 --threads 8 -x ${ref} -1 ${downdest}down_${cov}_C1776_S1_L001_R1_001.fastq -2 ${downdest}down_${cov}_C1776_S1_L001_R2_001.fastq -S ${aligndest}cneo_1776_${cov}.sam &> ${aligndest}1776_${cov}_log.txt &
bowtie2 --threads 8 -x ${ref} -1 ${downdest}down_${cov}_C1939_S2_L001_R1_001.fastq -2 ${downdest}down_${cov}_C1939_S2_L001_R2_001.fastq -S ${aligndest}cneo_1939_${cov}.sam &> ${aligndest}1939_${cov}_log.txt &
bowtie2 --threads 8 -x ${ref} -1 ${downdest}down_${cov}_C2898_S3_L001_R1_001.fastq -2 ${downdest}down_${cov}_C2898_S3_L001_R2_001.fastq -S ${aligndest}cneo_2898_${cov}.sam &> ${aligndest}2898_${cov}_log.txt

samtools view -bS ${aligndest}cneo_1776_${cov}.sam > ${aligndest}cneo_1776_${cov}.bam &
samtools view -bS ${aligndest}cneo_1939_${cov}.sam > ${aligndest}cneo_1939_${cov}.bam &
samtools view -bS ${aligndest}cneo_2898_${cov}.sam > ${aligndest}cneo_2898_${cov}.bam

samtools sort -o ${aligndest}cneo_1776_${cov}.sorted.bam ${aligndest}cneo_1776_${cov}.bam &
samtools sort -o ${aligndest}cneo_1939_${cov}.sorted.bam ${aligndest}cneo_1939_${cov}.bam &
samtools sort -o ${aligndest}cneo_2898_${cov}.sorted.bam ${aligndest}cneo_2898_${cov}.bam

echo Alignment done------------------------------





##Analysis of grubii variant.
refa=/mithril/Data/NGS/Reference/cneo/grubii/grubii.fa
piledest=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample${cov}X/mpileup/
mkdir ${piledest}

echo Consensus------------------------------------
##samtools faidx done on ref
samtools mpileup -uf $refa ${aligndest}cneo_1776_${cov}.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${piledest}cneo_1776_${cov}.cons.fq &
samtools mpileup -uf $refa ${aligndest}cneo_1939_${cov}.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${piledest}cneo_1939_${cov}.cons.fq &
samtools mpileup -uf $refa ${aligndest}cneo_2898_${cov}.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${piledest}cneo_2898_${cov}.cons.fq &


samtools mpileup -uf $refa ${aligndest}cneo_1776_${cov}.sorted.bam | bcftools call -c - | vcfutils.pl varFilter -D100 > ${piledest}cneo_1776_${cov}.vcf &
samtools mpileup -uf $refa ${aligndest}cneo_1939_${cov}.sorted.bam | bcftools call -c - | vcfutils.pl varFilter -D100 > ${piledest}cneo_1939_${cov}.vcf &
samtools mpileup -uf $refa ${aligndest}cneo_2898_${cov}.sorted.bam | bcftools call -c - | vcfutils.pl varFilter -D100 > ${piledest}cneo_2898_${cov}.vcf


##convert to fasta to appease parsnp
mkdir ${piledest}fastas
seqtk seq -a ${piledest}cneo_1776_${cov}.cons.fq > ${piledest}fastas/cneo_1776_${cov}.fasta
seqtk seq -a ${piledest}cneo_1939_${cov}.cons.fq > ${piledest}fastas/cneo_1939_${cov}.fasta
seqtk seq -a ${piledest}cneo_2898_${cov}.cons.fq > ${piledest}fastas/cneo_2898_${cov}.fasta
echo Consensus done-------------------------------






echo Parsnp---------------------------------------
#spades=/atium/Data/NGS/Aligned/160921_SSSSS_Crypt/grubii/downsample30X/spades/
parsnpdest=~/Dropbox/Lab/fungus_zhang/

~/software/parsnp/parsnp -r $refa -d ${piledest}/fastas -p 12 -o ${parsnpdest}parsnp_${cov}_cons
#~/software/parsnp/parsnp -r $refa -d $spades -p 12 -o ${out}parsnp_${cov}_spades

mv ${parsnpdest}parsnp_${cov}_cons/parsnp.ggr ${parsnpdest}parsnp_${cov}_cons/parsnp_cons_${cov}.ggr
#mv ${parsnpdest}parsnp_${cov}_spades/parsnp.ggr ${parsnpdest}parsnp_${cov}_spades/parsnp_spades.ggr
