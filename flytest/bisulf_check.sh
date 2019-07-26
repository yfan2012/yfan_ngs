#!/bin/bash

fqdir=~/data/190718_flytest/reads
datadir=~/data/190718_flytest
trimdir=~/data/190718_flytest/trimmed

if [ $1 == prep_genomes ] ; then
    ##for i in DRME DRME_2 human ecoli;
    for i in elegans ;
    do 
       bismark_genome_preparation --bowtie2 --genomic_composition --parallel 36 $datadir/ref/$i
    done
fi

if [ $1 == trim ] ; then
    for i in meth unmeth ;
    do
	mkdir -p $datadir/trimmed_$i
	trim_galore -q 30 --paired $fqdir/${i}*_R1_*.fastq.gz $fqdir/${i}*_R2_*.fastq.gz -o $datadir/trimmed_$i &> $datadir/trimmed_$i/trimlog.txt
    done
fi

if [ $1 == trimmomatic ] ; then
    ##try some v aggressive trimming
    for i in meth unmeth ;
    do
	mkdir -p $datadir/trimmomatic
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 -phred33 \
	     $fqdir/${i}*_R1_*.fastq.gz $fqdir/${i}*_R2_*.fastq.gz \
	     $datadir/trimmomatic/${i}_fwd_paired.fq.gz $datadir/trimmomatic/${i}_fwd_unpaired.fq.gz \
	     $datadir/trimmomatic/${i}_rev_paired.fq.gz $datadir/trimmomatic/${i}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:all.fa:4:10:10 LEADING:32 TRAILING:32 SLIDINGWINDOW:4:30 MINLEN:36 AVGQUAL:32
    done
fi

if [ $1 == btidx ] ; then
    for i in  DRME DRME_2 human ecoli ;
    do
	prefix=`basename $datadir/ref/$i/*.fa .fa`
	bowtie2-build --threads 36 $datadir/ref/$i/*.fa $datadir/ref/$i/$prefix
    done
fi


if [ $1 == bismark ] ; then
    for i in meth unmeth ;
    do
	##for genome in DRME DMRE_2 human ecoli ;
	for genome in elegans ;
	do
	    mkdir -p $datadir/${i}_$genome
	    bismark --bam --bowtie2 \
		    -p 36 \
		    --score_min L,0,-0.6 \
		    --un \
		    --non_directional \
		    --prefix ${i}_$genome \
		    $datadir/ref/$genome \
		    -1 ${trimdir}_$i/${i}*val_1.fq.gz \
		    -2 ${trimdir}_$i/${i}*val_2.fq.gz \
		    -o $datadir/${i}_$genome \
		&> $datadir/${i}_$genome/bismarklog.txt
	done
    done
fi


if [ $1 == bismark_r1 ] ; then
    for i in meth unmeth ;
    do
	trimmodir=$datadir/trimmomatic
	##for genome in DRME DMRE_2 human ecoli ;
	for genome in DRME_2 ;
	do
	    r1=`ls ${trimmodir}/${i}_fwd_paired.fq.gz`
	    r2=`ls ${trimmodir}/${i}_rev_paired.fq.gz`
	    mkdir -p $datadir/${i}_r1_$genome
	    bismark --bam --bowtie2 \
		    -p 36 \
		    --score_min L,0,-0.6 \
		    --un \
		    --non_directional \
		    --prefix ${i}_r1_$genome \
		    $datadir/ref/$genome \
		    --se $r1,$r2 \
		    -o $datadir/${i}_r1_$genome \
		&> $datadir/${i}_r1_$genome/bismarklog.txt
	done
    done
fi
		      
		      
		      
if [ $1 == meth_extract ] ; then
    for i in meth unmeth ;
    do
	bismark_methylation_extractor --parallel 36 -p --include_overlap -o $datadir/$i $datadir/$i/$i*.bam
    done
fi

    
if [ $1 == trim_nextera ] ; then
    ##try some v aggressive trimming of the gdna data
    fqdir=/kyber/Data/NGS/Raw/190725_flywtf
    datadir=~/data/190725_flywtf
    mkdir -p $datadir
    for i in meth unmeth ;
    do
	mkdir -p $datadir/trimmomatic
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 -phred33 \
	     $fqdir/${i}*_R1_*.fastq.gz $fqdir/${i}*_R2_*.fastq.gz \
	     $datadir/trimmomatic/${i}_fwd_paired.fq.gz $datadir/trimmomatic/${i}_fwd_unpaired.fq.gz \
	     $datadir/trimmomatic/${i}_rev_paired.fq.gz $datadir/trimmomatic/${i}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:all.fa:4:10:10 LEADING:32 TRAILING:32 SLIDINGWINDOW:4:30 MINLEN:36 AVGQUAL:32
    done
fi
