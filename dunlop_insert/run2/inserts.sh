#!/bin/bash

datadir=/mithril/Data/NGS/projects/dunlop_insert/run2
ref=$datadir/reference/construct.fa

if [ $1 == find ] ; then
    mkdir -p $datadir/positions

    for i in NT284 NT285 NT286 NT287 ;
    do
	python ~/Code/yfan_ngs/dunlop_insert/find_insert.py \
	       -b $datadir/align/$i.sorted.bam \
	       -o $datadir/positions/$i.positions.csv
    done
fi
