#!/bin/bash


##set up some direcotry names
assemblydir=/dilithium/Data/NGS/Aligned/171115_sears_fmt/assemblies/all
outdir=~/Dropbox/yfan/sears_fmt/ani_new
mkdir -p $outdir

##use the python package pyani to calculate percent identity between all samples
average_nucleotide_identity.py -v -f -o $outdir -i $assemblydir -g --gformat png,pdf,eps,svg --write_excel &> $outdir/ani_log.txt
