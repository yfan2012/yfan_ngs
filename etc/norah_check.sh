#!/bin/bash

datadir=~/work/ngs/190215_NA12878_cDNA
refdir=/home-4/yfan7@jhu.edu/work/Reference/human/hg38
if [ $1 == classify ] ; then
    ml gcc
    dbdir=~/scratch/centrifuge_db
    srcdir=~/scratch/centrifuge
    outdir=$datadir/classification
    
    mkdir -p $outdir

    prefix=NA12878
    
    $srcdir/centrifuge -p 36 -x $dbdir/abvm -1 $datadir/fastqs/*R1*fastq.gz -2 $datadir/fastqs/*R2*.fastq.gz -S $datadir/classification/$prefix.txt --report-file $datadir/classification/${prefix}_report.tsv
    $srcdir/centrifuge-kreport -x $dbdir/abvm $datadir/classification/$prefix.txt > $datadir/classification/kreport_$prefix.txt
fi


if [ $1 == align ] ; then
    ml samtools
    
    refdir=/home-4/yfan7@jhu.edu/work/Reference/human/hg38
    aligndir=$datadir/align

    mkdir -p $aligndir
    
    bowtie2 -p 36 -x $refdir/GRCH38 -1 $datadir/fastqs/*R1*fastq.gz -2 $datadir/fastqs/*R2*.fastq.gz | samtools view -bS - | samtools sort -o $aligndir/NA12878.sorted.bam
    samtools index $aligndir/NA12878.sorted.bam
fi

if [ $1 == align_unpaired ] ; then
    ml samtools
    

    aligndir=$datadir/align

    mkdir -p $aligndir
    
    bowtie2 -p 36 -x $refdir/GRCH38 -U $datadir/fastqs/NA12878-cDNA_all.fastq.gz | samtools view -bS - | samtools sort -o $aligndir/NA12878_unpaired.sorted.bam
    samtools index $aligndir/NA12878_unpaired.sorted.bam
fi


if [ $1 == trim ] ; then
    mkdir -p $datadir/trimmed
    
    java -jar ~/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 36 -phred33 $datadir/fastqs/*R1*.fastq.gz $datadir/fastqs/*R2*.fastq.gz $datadir/trimmed/NA12878_forward_paired.fq.gz $datadir/trimmed/NA12878_forward_unpaired.fq.gz $datadir/trimmed/NA12878_reverse_paired.fq.gz $datadir/trimmed/NA12878_reverse_unpaired.fq.gz ILLUMINACLIP:nebnextultra2.fa:6:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fi

if [ $1 == align_trimmed ] ; then
    ml samtools
    
    refdir=/home-4/yfan7@jhu.edu/work/Reference/human/hg38
    aligndir=$datadir/align

    mkdir -p $aligndir
    
    bowtie2 -p 36 -x $refdir/GRCH38 -1 $datadir/trimmed/*forward_paired*fq.gz -2 $datadir/trimmed/*reverse_paired*.fq.gz | samtools view -bS - | samtools sort -o $aligndir/NA12878_trimmed_paired.sorted.bam
    samtools index $aligndir/NA12878_trimmed_paired.sorted.bam
fi


if [ $1 == hisat_trimmed ] ; then
    ml samtools

    refdir=/home-4/yfan7@jhu.edu/work/Reference/hisat
    aligndir=$datadir/align
    ##hisat2-build $refdir/GRCH38.fa $refdir/GRCH38
    hisat2 -p 36 -q -x $refdir/GRCH38 -1 $datadir/trimmed/*forward_paired*fq.gz -2 $datadir/trimmed/*reverse_paired*.fq.gz | samtools view -bS - | samtools sort -o $aligndir/NA12878_hisat_trimmed_paired.sorted.bam
    samtools index $aligndir/NA12878_hisat_trimmed_paired.sorted.bam
fi

if [ $1 == hisat_trimmed_worm ] ; then
    ml samtools
    
    refdir=/home-4/yfan7@jhu.edu/work/Reference/celegans
    aligndir=$datadir/align
    
    hisat2-build $refdir/celegans.fa $refdir/celegans
    hisat2 -p 36 -q -x $refdir/celegans -1 $datadir/trimmed/*forward_paired*fq.gz -2 $datadir/trimmed/*reverse_paired*.fq.gz | samtools view -bS - | samtools sort -o $aligndir/NA12878_hisat_trimmed_paired_celegans.sorted.bam
    samtools index $aligndir/NA12878_hisat_trimmed_paired_celegans.sorted.bam
fi

if [ $1 == align_trimmed_worm ] ; then
    ml samtools
    
    refdir=/home-4/yfan7@jhu.edu/work/Reference/celegans
    aligndir=$datadir/align

    mkdir -p $aligndir

    bowtie2-build $refdir/celegans.fa $refdir/celegans
    bowtie2 -p 36 -x $refdir/celegans -1 $datadir/trimmed/*forward_paired*fq.gz -2 $datadir/trimmed/*reverse_paired*.fq.gz | samtools view -bS - | samtools sort -o $aligndir/NA12878_trimmed_paired_celegans.sorted.bam
    samtools index $aligndir/NA12878_trimmed_paired_celegans.sorted.bam
fi

if [ $1 == classify_trimmed ] ; then
    ml gcc
    dbdir=~/scratch/centrifuge_db
    srcdir=~/scratch/centrifuge
    outdir=$datadir/classification
    
    mkdir -p $outdir

    prefix=NA12878
    
    $srcdir/centrifuge -p 36 -x $dbdir/abvm -1 $datadir/trimmed/*forward_paired*fq.gz -2 $datadir/trimmed/*reverse_paired*.fq.gz -S $datadir/classification/$prefix.trimmed.txt --report-file $datadir/classification/${prefix}_report.trimmed.tsv
    $srcdir/centrifuge-kreport -x $dbdir/abvm $datadir/classification/$prefix.trimmed.txt > $datadir/classification/kreport_$prefix.trimmed.txt
fi
