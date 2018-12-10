#!/bin/bash



if [ $1 == old ] ; then
    assembledir=~/work/ngs/171115_sears_fmt/assemble
    
    alldir=~/work/ngs/171115_sears_fmt/assemblies
    mkdir -p $alldir
    
    for i in $assembledir/* ; do
	prefix=${i#$assembledir/}
	cp $i/contigs.fasta $alldir/$prefix.spades.fasta
    done
fi


if [ $1 == newrun ] ; then
    assembledir=~/work/ngs/180628_sears_fmt2/assemble

    alldir=~/work/ngs/171115_sears_fmt/assemblies
    mkdir -p $alldir

    for i in $assembledir/* ; do
	prefix=${i#$assembledir/}
	cp $i/contigs.fasta $alldir/${prefix}_new.spades.fasta
    done
fi
