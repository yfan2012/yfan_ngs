#!/bin/bash

rawdir=/mithril/Data/NGS/Raw/200918_dunlop_nathan
datadir=/mithril/Data/NGS/projects/dunlop_insert

ref=$datadir/reference/construct.fa

if [ $1 == find ] ; then
    mkdir -p $datadir/positions

    for i in B1_S1 B2_S2 B3_S3 B4_S4 ;
    do
	python ~/Code/yfan_ngs/dunlop_insert/find_insert.py \
	       -b $datadir/align/$i.sorted.bam \
	       -o $datadir/positions/$i.positions.csv
    done
fi
