#!/bin/bash

datadir=/kyber/Data/NGS/projects/190513_hardwick
outdir=~/Dropbox/yfan/hardwick

if [ $1 == maker_ctlfiles ] ; then
    maker -CTL
fi


if [ $1 == augustus ] ; then
    ##thankfully, augustus has crypto already built in
    
    augdir=$datadir/augustus
    mkdir -p $augdir

    for i in 1694 178 197 6341 ;
    do
	augustus --species=cryptococcus $datadir/longtigs/${i}_over20k.fasta > $augdir/${i}_over20k.gff &
	augustus --species=cryptococcus_neoformans_gattii $datadir/longtigs/${i}_over20k.fasta > $augdir/${i}_over20k_gattii.gff &
	augustus --species=cryptococcus_neoformans_neoformans_B $datadir/longtigs/${i}_over20k.fasta > $augdir/${i}_over20k_B.gff &
	augustus --species=cryptococcus_neoformans_neoformans_JEC21 $datadir/longtigs/${i}_over20k.fasta > $augdir/${i}_over20k_JEC21.gff &
    done
fi
    
if [ $1 == parsnp ] ; then
    ##set each patient isolate as reference
    for i in 1694 178 197 6341 ;
    do
	mkdir -p $outdir/parsnp_$i
	parsnp -p 12 -r $outdir/assemblies/$i.fasta -d $outdir/assemblies -o $outdir/parsnp_$i
	harvesttools -i $outdir/parsnp_$i/parsnp.ggr -V $outdir/parsnp_$i/strain_snps.vcf
    done
fi

if [ $1 == countsnps ] ; then
    for i in 1694 178 197 6341 ;
    do
	Rscript ~/Code/utils/count_snps.R -i $outdir/parsnp_$i/strain_snps.vcf -o $outdir/parsnp_$i
    done
fi
