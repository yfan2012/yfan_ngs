#!/bin/bash


##set up some direcotry names
assemblydir=/dilithium/Data/NGS/Aligned/171115_sears_fmt/assemblies/all
outdir=~/Dropbox/Lab/sears_fmt/ani


##use the python package pyani to calculate percent identity between all samples
average_nucleotide_identity.py -v -f -o $outdir -i $assemblydir -g --gformat png,pdf,eps,svg --write_excel &> $outdir/ani_log.txt
