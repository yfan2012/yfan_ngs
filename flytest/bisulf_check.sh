#!/bin/bash

datadir=/kyber/Data/NGS/projects/190718_flytest
fqdir=/kyber/Data/NGS/Raw/190718_flytest
trimdir=$datadir/trimmomatic

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
    cat ~/software/Trimmomatic-0.39/adapters/*.fa > ~/Code/yfan_ngs/flytest/all.fa
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
	for genome in DRME ;
	do
	    r1=`ls ${trimmodir}/${i}_fwd_paired.fq.gz`
	    r2=`ls ${trimmodir}/${i}_rev_paired.fq.gz`
	    bismark --bam --bowtie2 \
		    -p 36 \
		    --score_min L,0,-0.6 \
		    --un \
		    --non_directional \
		    --prefix ${i}_se_$genome \
		    $datadir/ref/$genome \
		    --se $r1,$r2 \
		    -o $datadir/$i \
		&> $datadir/$i/bismarklog.txt
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

if [ $1 == centrifuge ] ; then
    ##run cent on bigmem, obvs
    ml gcc
    centdir=~/scratch/centrifuge
    dbdir=~/scratch/centrifuge_db
    datadir=~/work/ngs/190725_flywtf

    mkdir -p $datadir/classification

    $centdir/centrifuge -p 36 -x $dbdir/abvm \
			-1 $datadir/meth_fwd_paired.fq.gz,$datadir/unmeth_fwd_paired.fq.gz \
			-2 $datadir/meth_rev_paired.fq.gz,$datadir/unmeth_rev_paired.fq.gz \
			-S $datadir/classification/all.txt --report-file $datadir/classification/all_report.tsv
    $centdir/centrifuge-kreport -x $dbdir/abvm $datadir/classification/all.txt > $datadir/classification/all_kreport.txt
fi

if [ $1 == btidx ] ; then
    datadir=~/data/190718_flytest
    refs=$datadir/ref
    
    for i in DRME_A7 elegans ;
    do
	bowtie2-build --threads 36 $datadir/ref/$i/$i.fa $datadir/ref/$i/$i
    done
fi    

if [ $1 == align_errthing ] ; then
    refs=$datadir/ref
    new=~/data/190725_flywtf
    newfq=$new/trimmomatic
    aligndir=$new/align
    mkdir -p $aligndir
    
    for i in $refs/*/*.fa ;
    do
	prefix=`basename $i .fa`
	bowtie2 -p 36 -x $refs/$prefix/$prefix \
		-1 $newfq/meth_fwd_paired.fq.gz \
		-2 $newfq/meth_rev_paired.fq.gz \
	    | samtools view -bS - | samtools sort -o $aligndir/$prefix.meth.sorted.bam
	samtools index $aligndir/$prefix.meth.sorted.bam
    done
fi
    

if [ $1 == makeblastdb ] ; then
    refdir=$datadir/ref/DRME
    makeblastdb -in $refdir/DRME.fa -out $refdir/DRMEdb -dbtype nucl
fi

    
if [ $1 == blast ] ; then
    readdir=~/data/190725_flywtf/trimmomatic
    refdir=$datadir/ref/DRME
    seqtk seq -a $readdir/meth_fwd_paired.fq.gz > $readdir/meth_fwd_paired.fa

    mkdir -p ~/data/190725_flywtf/blast
    blastn -query $readdir/meth_fwd_paired.fa -db $refdir/DRMEdb -outfmt 7 -num_threads 48 -out ~/data/190725_flywtf/blast/DRME.tsv
fi
