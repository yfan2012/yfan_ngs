#!/bin/bash

index=SI-TT-A10
samp=human_skin
run=220307_NrlSkin
datadir=/dilithium/Data/NGS/projects/human_skin_yfan

dbdir=/dilithium/Data/NGS/Reference/kraken2


##Try to set up kraken nt database
if [ $1 == dldb ] ; then
    ##the dustmasker step takes like two days
    kraken2-build \
	--download-library nt \
	--threads 36 \
	--use-ftp \
	--db $dbdir/nt
fi

if [ $1 == dltax ] ; then
    kraken2-build \
	--download-taxonomy \
	--db $dbdir/nt
fi

if [ $1 == builddb ] ; then
    kraken2-build \
	--build \
	--threads 36 \
	--db $dbdir/nt
fi


##cellranger pipeline - not sure what the point of mkfastq is - let's try with and without
if [ $1 == sampsheet ] ; then
    ##make sample sheet
    touch $datadir/raw/$run/$run.csv
    echo Lane,Sample,Index >> $datadir/raw/$run/$run.csv
    echo 1,human_skin,SI-TT-A10 >> $datadir/raw/$run/$run.csv
fi


if [ $1 == mkfastq ] ; then
    cellranger mkfastq \
	       --id=mkfastq \
	       --sample-sheet=$datadir/raw/$run/$run.csv \
	       --run=$datadir/raw/$run/Files
    mv ./mkfastq/ $datadir/
fi   

ref=/dilithium/Data/NGS/Reference/cellranger/refdata-gex-GRCh38-2020-A
if [ $1 == count ] ; then
    cellranger count \
	       --id=count \
	       --fastqs=$datadir/mkfastq/outs/fastq_path \
	       --transcriptome=$ref
    mv ./count/ $datadir/
fi
