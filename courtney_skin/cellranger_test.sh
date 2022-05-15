#!/bin/bash

samp=human_skin
run=cellranger-tiny-bcl-1.2.0
datadir=/dilithium/Data/NGS/projects/human_skin_yfan/test


if [ $1 == mkfastq ] ; then
    cellranger mkfastq \
	       --id=mkfastq \
	       --output-dir=$datadir/bcl \
	       --sample-sheet=$datadir/raw/cellranger-tiny-bcl-simple-1.2.0.csv \
	       --run=$datadir/raw/$run
    mv ./mkfastq/ $datadir/
fi   


if [ $1 == count ] ; then
    cellranger count \
	       --id=count \
	       --fastqs=$datadir/raw/pbmc_1k_v3_fastqs \
	       --transcriptome=$datadir/Reference/refdata-cellranger-GRCh38-3.0.0
    mv ./count/ $datadir/
fi
