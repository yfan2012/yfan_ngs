#!/bin/bash

projdir=/mithril/Data/NGS/projects/dunlop_insert
datadir=$projdir/run3
ref1=$projdir/refs/construct1.fa
ref2=$projdir/refs/construct2.fa

if [ $1 == test ] ; then
    mkdir -p $datadir/exact
    i=NT278
    python ~/Code/yfan_ngs/dunlop_insert/find_exact.py \
	   -1 $datadir/trimmed/${i}_fwd_paired.fq.gz \
	   -2 $datadir/trimmed/${i}_rev_paired.fq.gz \
	   -3 $datadir/trimmed/${i}_fwd_unpaired.fq.gz \
	   -4 $datadir/trimmed/${i}_rev_unpaired.fq.gz \
	   -r $projdir/refs/construct1.fa \
	   -o $datadir/exact/test.positions.csv \
	   -t 4
fi

if [ $1 == findmag ] ; then
    mkdir -p $datadir/exact

    for i in NT278 NT296 NT297 ;
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
fi
if [ $1 == findilid ] ; then
    for i in NT279 NT298 NT299 ;
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
