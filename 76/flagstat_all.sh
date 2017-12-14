#!/bin/bash

datadir=/dilithium/Data/NGS/Aligned/170928_fungus76
mkdir -p $datadir/stat_grubii
mkdir -p $datadir/stat_jec21
mkdir -p $datadir/stat_b-3501a

for samp in $datadir/align_grubii/*.sorted.bam ;
do
    name=` echo ${samp%.sorted.bam} | cut -d '/' -f 8 `
    echo $name
    samtools flagstat $samp > $datadir/stat_grubii/$name.stat.txt &
    samtools flagstat $datadir/align_jec21/$name.sorted.bam > $datadir/stat_jec21/$name.stat.txt &
    samtools flagstat $datadir/align_b-3501a/$name.sorted.bam > $datadir/stat_b-3501a/$name.stat.txt
done
