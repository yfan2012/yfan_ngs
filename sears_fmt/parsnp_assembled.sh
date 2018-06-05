#!/bin/bash

datadir=/dilithium/Data/NGS/Aligned/171115_sears_fmt
assembledir=$datadir/assemblies
ref=/mithril/Data/NGS/Reference/ecoli/ecoli.fasta


outdir=~/Dropbox/yfan/sears_fmt/parsnp_assembly
mkdir -p $outdir

if [ $1 == ecoli ] ; then
    ##parsnp -r $ref -d $assembledir -p 12 -o $outdir -c
    
    ##ecoli separated from kleb by hand
    ecolidir=$assembledir/ecoli
    klebdir=$assembledir/kleb
    
    ecoliout=~/Dropbox/Lab/sears_fmt/parsnp_ecoli
    mkdir -p $ecoliout
    ##parsnp -r ! -d $ecolidir -p 12 -o $ecoliout -c
    
    
    
    klebref=/mithril/Data/NGS/Reference/kpneumo/NC_016845.fasta
    klebout=~/Dropbox/Lab/sears_fmt/parsnp_kleb
    mkdir -p $klebout
    parsnp -r ! -d $klebdir -p 12 -o $klebout -c
fi



if [ $1 == ecolidonor ] ; then
    ##make a different tree for each donor, comparing all recipients to a single donor
    recipdir=$assembledir/ecoli/recipient
    outdir=~/Dropbox/Lab/sears_fmt

    for i in $assembledir/ecoli/donor/*fasta ;
    do
	dsamp=`basename $i .spades.fasta`
	doutdir=$outdir/parsnp_$dsamp
	mkdir -p $doutdir
	parsnp -r $i -d $recipdir -p 12 -o $doutdir -c
    done
fi

if [ $1 == ecolidonor2 ] ; then
    ##make a different tree for each donor, comparing all recipients to a single donor
    ##Eject 1R12C and 3R12C    
    recipdir=$assembledir/ecoli/recipient2
    outdir=~/Dropbox/Lab/sears_fmt

    for i in $assembledir/ecoli/donor/*fasta ;
    do
	dsamp=`basename $i .spades.fasta`
	doutdir=$outdir/parsnp_${dsamp}2
	mkdir -p $doutdir
	parsnp -r $i -d $recipdir -p 12 -o $doutdir -c
    done
fi
       


if [ $1 == all ] ; then
    #try parsnp for all bacteria
    outdir=~/Dropbox/yfan/sears_fmt/parsnp_all
    alldir=$datadir/assemblies/all
    mkdir -p outdir
    
    parsnp -r $alldir/ecoli_NC101.fa -d $alldir -p 12 -o $outdir -c
fi


if [ $1 == withref ] ; then
	outdir=~/Dropbox/yfan/sears_fmt/parsnp_ecoli
	gendir=$datadir/assemblies/ecoli/all
	mkdir -p $outdir
	parsnp -r $gendir/ecoli_NC101.fa -d $gendir -p 12 -o $outdir -c
	harvesttools -i $outdir/parsnp.ggr -V $outdir/ecoli.vcf
	
fi



