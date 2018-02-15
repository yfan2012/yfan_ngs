#!/bin/bash

##Run centrifuge all samps to check they're all ecoli. I think at least 4 are not.

##use bacteria/virus/archea database
dbdir=~/scratch/centrifuge_db
srcdir=~/scratch/centrifuge

datadir=~/work/ngs/171115_sears_fmt/fastqs
outdir=~/work/ngs/171115_sears_fmt/classification

mkdir -p $outdir

for i in $datadir/*R1_001.fastq.gz ;
do
    read2=${i%R1_001.fastq.gz}R2_001.fastq.gz

    file=${i#$datadir/}
    prefix=${file%_L001_R1_001.fastq.gz}

    mkdir -p $datadir/classification/$prefix
    
    $srcdir/centrifuge -p 36 -x $dbdir/abv -1 $i -2 $read2 -S $datadir/classification/$prefix.txt --report-file $datadir/classification/${prefix}_report.tsv
    $srcdir/centrifuge-kreport -x $dbdir/abv $datadir/classification/$prefix.txt > $datadir/classification/kreport_$prefix.txt
done
