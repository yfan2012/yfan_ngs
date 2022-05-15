#!/bin/bash

samp=SI-TT-A10
projdir=/dilithium/Data/NGS/projects/human_skin
rawdir=$projdir/Samples/$samp/Files
datadir=${projdir}_yfan

dbdir=/dilithium/Data/NGS/Reference/kraken2

if [ $1 == dldb ] ; then
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
