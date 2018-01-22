#!/bin/bash

assembledir=~/work/ngs/171115_sears_fmt/assemble

alldir=~/work/ngs/171115_sears_fmt/assemblies
mkdir -p $alldir

for i in $assembledir/* ; do
    prefix=${i#$assembledir/}
    cp $i/contigs.fasta $alldir/$prefix.spades.fasta
done
