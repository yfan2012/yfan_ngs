#!/bin/bash

datadir=/atium/Data/projects/ambic_cho_yfan

if [ $1 == getref ] ; then
    ##refs from ncbi
    mkdir -p $datadir/ref
    
    wget -O $datadir/ref/CHO.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/135/GCF_000223135.1_CriGri_1.0/GCF_000223135.1_CriGri_1.0_genomic.fna.gz

    wget -O $datadir/ref/CHO.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/135/GCF_000223135.1_CriGri_1.0/GCF_000223135.1_CriGri_1.0_genomic.gff.gz

    
    gunzip -c $datadir/ref/CHO.gff.gz > $datadir/ref/CHO.gff
    gffread $datadir/ref/CHO.gff -T -o $datadir/ref/CHO.gtf

    ##unzip fa for fun
    gunzip -c $datadir/ref/CHO.fa.gz > $datadir/ref/CHO.fa
fi


if [ $1 == addigg ] ; then
    cp /pym/Data/reference/cho/ambic_sigma_IgG.fa $datadir/ref/
    
    cat $datadir/ref/ambic_sigma_IgG.fa >> $datadir/ref/CHO.fa

    ##gtf manually written from kevin (find in slack dms from alice)
    ##need to make sure everything in there is tabs not spaces - did this manually
    cat $datadir/ref/ambic_sigma_IgG.gtf >> $datadir/ref/CHO.gtf
fi

if [ $1 == getgtf ] ; then
    ##gffread conversion to gtf doesn't preserve gene ids which cell ranger complains about
    ##try downloading from ncbi to see if it's helpful
    
    wget -O $datadir/ref/CHO.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/223/135/GCF_000223135.1_CriGri_1.0/GCF_000223135.1_CriGri_1.0_genomic.gtf.gz
    gunzip -c $datadir/ref/CHO.gtf.gz > $datadir/ref/CHO.gtf
fi

if [ $1 == gtfigg ] ; then
    cat $datadir/ref/ambic_sigma_IgG.gtf >> $datadir/ref/CHO.gtf
fi


if [ $1 == liftoff ] ; then
    ##see if we can transfer good annotation to the picr genome so we know what genes they are
    mkdir -p $datadir/liftoff
    cp /atium/Data/projects/ambic_cho/cho_ref/CriGri.PICRH_genome_igg/picr.CriGri-PICR_igg/fasta/genome.fa $datadir/liftoff/genome.fa

    liftoff \
	-g $datadir/ref/CHO_annot.gff \
	-o $datadir/liftoff/CHO_picr_liftover.gff \
	-u $datadir/liftoff/unmapped_featrues.txt \
	-infer_genes \
	$datadir/liftoff/genome.fa \
	$datadir/ref/CHO.fa
fi


if [ $1 == liftoffgtf ] ; then
    gffread $datadir/liftoff/CHO_picr_liftover.gff -T -o $datadir/liftoff/CHO_picr_liftover.gtf
fi

