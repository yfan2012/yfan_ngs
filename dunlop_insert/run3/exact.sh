#!/bin/bash

projdir=/mithril/Data/NGS/projects/dunlop_insert
datadir=$projdir/run3

if [ $1 == test ] ; then
    mkdir -p $datadir/exact
    i=NT278
    python ~/Code/yfan_ngs/dunlop_insert/find_exact.py \
	   -1 $datadir/trimmed/${i}_fwd_paired.fq.gz \
	   -2 $datadir/trimmed/${i}_rev_paired.fq.gz \
	   -3 $datadir/trimmed/${i}_fwd_unpaired.fq.gz \
	   -4 $datadir/trimmed/${i}_rev_unpaired.fq.gz \
	   -r $projdir/refs/construct1.fa \
	   -o $datadir/exact/$i.positions.csv \
	   -t 4
fi
