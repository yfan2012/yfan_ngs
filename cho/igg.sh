#!/bin/bash

alicedir=/atium/Data/projects/ambic_cho
datadir=/atium/Data/projects/ambic_cho_yfan

if [ $1 == grabraw90 ] ; then
    mkdir -p $datadir/matrices
    
    rawdir=$alicedir/211104_first_sc_cho/schen239_176020/CHO_day90/outs/filtered_feature_bc_matrix
    cp $rawdir/*gz $datadir/matrices/
    
    for i in $datadir/matrices/*.gz ;
    do
	name=`basename $i .gz`
	mv $i $datadir/matrices/day90.$name.gz
    done
fi

if [ $1 == mkgtf ] ; then
    ##empty gene_id fields are all trna or rrna, so get rid of them
    grep -v "gene_id \"\"" $datadir/ref/CHO.gtf > $datadir/ref/CHO_mrna.gtf
    
    cellranger mkgtf \
	       $datadir/ref/CHO_mrna.gtf \
	       $datadir/ref/CHO.filtered.gtf \
               --attribute=gene_biotype:protein_coding
fi



if [ $1 == mkref ] ; then
    cellranger mkref \
	       --genome=CHO \
	       --fasta=$datadir/ref/CHO.fa \
	       --genes=$datadir/ref/CHO.filtered.gtf \
	       --nthreads=18 \
	       --memgb=128

    mv ./CHO $datadir/ref/
fi



if [ $1 == count0 ] ; then
    mkdir -p $datadir/cellranger
    cellranger count \
	       --id=CHO_day0_count \
	       --sample=CHO_day0 \
	       --transcriptome=$datadir/ref/CHO \
	       --fastqs=$alicedir/211104_first_sc_cho/schen239_176020/HN2WWDRXY/outs/fastq_path/HN2WWDRXY
			      
    mv ./CHO_day0_count $datadir/cellranger
fi

if [ $1 == count90 ] ; then
    mkdir -p $datadir/cellranger
    cellranger count \
	       --id=CHO_day90_count \
	       --sample=CHO_day90 \
	       --transcriptome=$datadir/ref/CHO \
	       --fastqs=$alicedir/211104_first_sc_cho/schen239_176020/HN2WWDRXY/outs/fastq_path/HN2WWDRXY
			      
    mv ./CHO_day90_count $datadir/cellranger
fi

##okay, this was fun and informative,
##but I'm switiching to alice's cellranger run for now
##to avoid any potential bugs in analysis
