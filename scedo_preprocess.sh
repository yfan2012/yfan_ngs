#!/bin/bash

samp=170412_scedo_hyphae
rawdir=/atium/Data/Nanopore/oxford/${samp}
dir=/atium/Data/Nanopore/Analysis/${samp}
rawtmp=/atium/Data/tmp/tmp_${samp}/


##explode tarball into rawtmp
if [ ! -d $rawtmp ]; then
    mkdir -p $rawtmp
    pigz -k -d $rawdir/${samp}_raw.tar.gz | tar -xf -C $rawtmp
fi
    
mkdir -p $dir

read_fast5_basecaller.py -i $rawtmp -t 10 -s $rawdir/${samp}_called/ -c FLO-MIN106_LSK108_linear.cfg --recursive &> $rawdir/albacore_log.txt 
    
tar --remove-files -cf $rawdir/${samp}_called.tar $rawdir/${samp}_called &&

pigz --best $rawdir/${samp}_called.tar



