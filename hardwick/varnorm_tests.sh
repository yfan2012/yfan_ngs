#!/bin/bash

datadir=/kyber/Data/NGS/projects/190513_hardwick/
allsamps='178 197 1694 6341'
ref=$datadir/reference/CRNE_H99.fasta


if [ $1 == trysamtools ] ; then
    for i in $allsamps ;
    do
	bcftools mpileup --threads 12 -Ou -f $ref $datadir/align/$i.sorted.bam \
	    | bcftools call --threads 12 -mv -Ov --ploidy 1 -o $datadir/varnorm/$i.bcftools.vcf
    done
fi

	
if [ $1 == tryvt ] ; then
    for i in $allsamps ;
    do
	vt normalize -r $ref -o $datadir/varnorm/$i.norm.vcf $datadir/varnorm/$i.bcftools.vcf
    done
fi

