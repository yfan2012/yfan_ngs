#!/bin/bash

rawdir=/kyber/Data/NGS/Raw/190513_hardwick
datadir=/kyber/Data/NGS/projects/190513_hardwick
outdir=~/Dropbox/yfan/hardwick

if [ $1 == gatherillumina ] ; then
    mkdir -p $rawdir
    for i in 1694 178 197 6341 ;
    do
	(cat ~/BaseSpace/Projects/190513_hardwick/Samples/$i/Files/$i*R1*.fastq.gz ~/BaseSpace/Projects/*190512_M01556_0154_000000000-C9C*/Samples/$i/Files/$i*R1*.fastq.gz > $rawdir/${i}_R1.fastq.gz
	 cat ~/BaseSpace/Projects/190513_hardwick/Samples/$i/Files/$i*R2*.fastq.gz ~/BaseSpace/Projects/*190512_M01556_0154_000000000-C9C*/Samples/$i/Files/$i*R2*.fastq.gz > $rawdir/${i}_R2.fastq.gz) &
    done
fi

if [ $1 == trim ] ; then
    mkdir -p $datadir
    mkdir -p $datadir/trimmed
    cp ~/software/Trimmomatic-0.39/adapters/NexteraPE-PE.fa ~/Code/yfan_ngs/hardwick/
    for i in 1694 178 197 6341 ;
    do
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 -phred33 \
	     $rawdir/${i}_R1.fastq.gz $rawdir/${i}_R2.fastq.gz \
	     $datadir/trimmed/${i}_fwd_paired.fq.gz $datadir/trimmed/${i}_fwd_unpaired.fq.gz \
	     $datadir/trimmed/${i}_rev_paired.fq.gz $datadir/trimmed/${i}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi

if [ $1 == spades ] ; then
    mkdir -p $datadir/spades
    for i in 1694 178 197 6341 ;
    do
	mkdir -p $datadir/spades/$i
	spades.py -1 $datadir/trimmed/${i}_fwd_paired.fq.gz -2 $datadir/trimmed/${i}_rev_paired.fq.gz -t 72 -m 300 -o $datadir/spades/$i
    done
fi

if [ $1 == nametigs ] ; then
    for i in 1694 178 197 6341 ;
    do
	mv $datadir/spades/$i/contigs.fasta $datadir/spades/$i/$i.fasta
    done
fi

if [ $1 == parsnp ] ; then
    mkdir -p $outdir/parsnp
    mkdir -p $outdir/assemblies
    for i in 1694 178 197 6341 ;
    do
	cp $datadir/spades/$i/$i.fasta $outdir/assemblies/
    done

    parsnp -p 12 -r ! -d $outdir/assemblies -o $outdir/parsnp
    harvesttools -i $outdir/parsnp/parsnp.ggr -V $outdir/parsnp/strain_snps.vcf
fi



if [ $1 == countsnps ] ; then
    Rscript ~/Code/utils/count_snps.R -i $outdir/parsnp/strain_snps.vcf -o $outdir/parsnp
fi



if [ $1 == dustmask ] ; then
    ##See if dustmasker can get rid of those 137 potentially mis-compared snps
    mkdir -p $datadir/dustmask
    for i in 1694 178 197 6341 ;
    do
	dustmasker -in $datadir/spades/$i/$i.fasta -out $datadir/dustmask/${i}_dustmasked.fasta -outfmt fasta &
    done
fi    
    
if [ $1 == parsnp_dm ] ; then
    mkdir -p $outdir/parsnp_dm
    mkdir -p $outdir/assemblies_dm
    for i in 1694 178 197 6341 ;
    do
	cp $datadir/dustmask/${i}_dustmasked.fasta $outdir/assemblies_dm/
    done

    parsnp -p 12 -r ! -d $outdir/assemblies_dm -o $outdir/parsnp_dm
    harvesttools -i $outdir/parsnp_dm/parsnp.ggr -V $outdir/parsnp_dm/strain_dm_snps.vcf
fi


if [ $1 == countsnps_dm ] ; then
    Rscript ~/Code/utils/count_snps.R -i $outdir/parsnp_dm/strain_dm_snps.vcf -o $outdir/parsnp_dm
fi

if [ $1 == nucmer ] ; then
    ##nucmer between the two patient strains

    mkdir -p $outdir/mummer
    nucmer -p $outdir/mummer/pt_strains $outdir/assemblies_dm/197_dustmasked.fasta $outdir/assemblies_dm/178_dustmasked.fasta

    mummerplot --filter --fat --png -p $outdir/mummer/pt_strains $outdir/mummer/pt_strains.delta

    dnadiff -p $outdir/mummer/pt_strains $outdir/assemblies_dm/197_dustmasked.fasta $outdir/assemblies_dm/178_dustmasked.fasta
fi
    
