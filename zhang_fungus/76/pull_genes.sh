#!/bin/bash

##Use gene references to pull out homologous seqs in the assemblies
##Organize them in a sensical way

datadir=/dilithium/Data/NGS/Aligned/170928_fungus76
refgenes=$datadir/ref_genes
mlstdir=$datadir/mlst
assemblydir=$datadir/assemble/parsnp_assemblies


if [ $1 == makedb ] ; then
    for i in $assemblydir/*fasta ; do
	samp=`basename $i .fasta`
	makeblastdb -in $i -out $assemblydir/$samp.db -dbtype nucl
    done
fi


if [ $1 == find ] ; then

    for ref in $refgenes/*.fasta ; do

	gene=`basename $ref .fasta`
	mkdir -p $mlstdir/$gene/regions
	
	for i in $assemblydir/*fasta ; do
	    samp=`basename $i .fasta`
	    blastn -query $ref -db $assemblydir/$samp.db -outfmt 7 -out $mlstdir/$gene/regions/$samp.tsv
	done
	
    done
fi
