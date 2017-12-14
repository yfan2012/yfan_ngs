#!/bin/bash

datadir=/dilithium/Data/NGS/Aligned/170928_fungus76

for i in $datadir/mpileup* ;
do
    for samp in $i/*.cons.fq ;
    do
	fasta=${samp%.cons.fq}.fasta
	seqtk seq -a $samp > $fasta
    done
done

