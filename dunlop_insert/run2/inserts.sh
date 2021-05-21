#!/bin/bash

projdir=/mithril/Data/NGS/projects/dunlop_insert
datadir=$projdir/run2
ref1=$projdir/refs/construct1.fa
ref2=$projdir/refs/construct2.fa

if [ $1 == test ] ; then
    mkdir -p $datadir/positions

    i=NT284
    python ~/Code/yfan_ngs/dunlop_insert/find_insert.py \
	   -b $datadir/align/$i.sorted.bam \
	   -r $ref1 \
	   -o $datadir/positions/$i.positions.csv
fi

if [ $1 == find_verobse ] ; then
    mkdir -p $datadir/positions

    for i in NT284 NT285 ;
    do
	python ~/Code/yfan_ngs/dunlop_insert/find_insert.py \
	       -b $datadir/align/$i.sorted.bam \
	       -r $ref1 \
	       -o $datadir/positions/$i.positions_verbose.csv \
	       -v
    done
    for i in NT286 NT287 ;
    do
	python ~/Code/yfan_ngs/dunlop_insert/find_insert.py \
	       -b $datadir/align/$i.sorted.bam \
	       -r $ref2 \
	       -o $datadir/positions/$i.positions_verbose.csv \
	       -v
    done
	
fi

if [ $1 == find ] ; then
    mkdir -p $datadir/positions

    for i in NT284 NT285 ;
    do
	python ~/Code/yfan_ngs/dunlop_insert/find_insert.py \
	       -b $datadir/align/$i.sorted.bam \
	       -r $ref1 \
	       -o $datadir/positions/$i.positions.csv
    done
    for i in NT286 NT287 ;
    do
	python ~/Code/yfan_ngs/dunlop_insert/find_insert.py \
	       -b $datadir/align/$i.sorted.bam \
	       -r $ref2 \
	       -o $datadir/positions/$i.positions.csv 
    done
	
fi
