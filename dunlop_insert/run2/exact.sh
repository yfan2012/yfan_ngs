#!/bin/bash

projdir=/mithril/Data/NGS/projects/dunlop_insert
datadir=$projdir/run2
ref1=$projdir/refs/construct1.fa
ref2=$projdir/refs/construct2.fa


if [ $1 == findmag ] ; then
    mkdir -p $datadir/exact

    for i in NT284 NT285 ;
    do
	(echo $i
	python ~/Code/yfan_ngs/dunlop_insert/find_exact.py \
           -1 $datadir/trimmed/${i}_fwd_paired.fq.gz \
           -2 $datadir/trimmed/${i}_rev_paired.fq.gz \
           -3 $datadir/trimmed/${i}_fwd_unpaired.fq.gz \
           -4 $datadir/trimmed/${i}_rev_unpaired.fq.gz \
           -r $ref1 \
           -o $datadir/exact/$i.positions.csv \
           -t 4) &
    done

    for i in NT286 NT287 ;
    do
	(echo $i
	python ~/Code/yfan_ngs/dunlop_insert/find_exact.py \
           -1 $datadir/trimmed/${i}_fwd_paired.fq.gz \
           -2 $datadir/trimmed/${i}_rev_paired.fq.gz \
           -3 $datadir/trimmed/${i}_fwd_unpaired.fq.gz \
           -4 $datadir/trimmed/${i}_rev_unpaired.fq.gz \
           -r $ref2 \
           -o $datadir/exact/$i.positions.csv \
           -t 4) &
    done

fi



if [ $1 == scarmag ] ; then
    mkdir -p $datadir/exact

    for i in NT284 NT285 ;
    do
	(echo $i
	python ~/Code/yfan_ngs/dunlop_insert/find_exact_scars.py \
           -1 $datadir/trimmed/${i}_fwd_paired.fq.gz \
           -2 $datadir/trimmed/${i}_rev_paired.fq.gz \
           -3 $datadir/trimmed/${i}_fwd_unpaired.fq.gz \
           -4 $datadir/trimmed/${i}_rev_unpaired.fq.gz \
           -r $ref1 \
           -o $datadir/exact/$i.scars.csv \
           -t 4) &
    done

    for i in NT286 NT287 ;
    do
	(echo $i
	python ~/Code/yfan_ngs/dunlop_insert/find_exact_scars.py \
           -1 $datadir/trimmed/${i}_fwd_paired.fq.gz \
           -2 $datadir/trimmed/${i}_rev_paired.fq.gz \
           -3 $datadir/trimmed/${i}_fwd_unpaired.fq.gz \
           -4 $datadir/trimmed/${i}_rev_unpaired.fq.gz \
           -r $ref2 \
           -o $datadir/exact/$i.scars.csv \
           -t 4) &
    done

fi
