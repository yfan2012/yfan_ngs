
#!/bin/bash

projdir=/mithril/Data/NGS/projects/dunlop_insert
datadir=$projdir/run3

if [ $1 == getref ] ; then
    ##dl ecoli reference
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -O $projdir/refs/ecoli.fa.gz
fi

if [ $1 == mergerefs ] ; then
    cat $projdir/refs/ptkei-dest-cre-loxp-sfgfp.fasta \
	| tr '[:lower:]' '[:upper:]' \
	| gzip > $projdir/refs/ptkei-dest-cre-loxp-sfgfp.fasta.gz
    
    cat $projdir/refs/ecoli.fa.gz $projdir/refs/ptkei-dest-cre-loxp-sfgfp.fasta.gz > $projdir/refs/ecoli_withplas.fasta.gz
fi

if [ $1 == mergeconstructs ] ; then
    cat $projdir/refs/construct1.fa $projdir/refs/ptkei-dest-cre-loxp-sfgfp.fasta \
	| tr '[:lower:]' '[:upper:]' \
	| gzip > $projdir/refs/construct1_plas.fa.gz

    cat $projdir/refs/construct2.fa $projdir/refs/ptkei-dest-cre-loxp-sfgfp.fasta \
	| tr '[:lower:]' '[:upper:]' \
	| gzip > $projdir/refs/construct2_plas.fa.gz
fi


if [ $1 == align ] ; then
    mkdir -p $datadir/check
    bwa index $projdir/refs/ecoli.fa.gz
    
    for i in NT278 NT279 NT296 NT297 NT298 NT299 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/ecoli.fa.gz \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_bwa_ecoli.sorted.bam
	samtools index $datadir/check/${i}_bwa_ecoli.sorted.bam
    done
fi

if [ $1 == bowtie ] ; then
    bowtie2-build $projdir/refs/ecoli.fa.gz $projdir/refs/ecoli
    for i in NT278 NT279 NT296 NT297 NT298 NT299 ;
    do
	bowtie2 \
	    -p 36 \
	    -x $projdir/refs/ecoli \
	    -1 $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    -2 $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_bt2_ecoli.sorted.bam
	samtools index $datadir/check/${i}_bt2_ecoli.sorted.bam
    done
fi

if [ $1 == align_withplas ] ; then
    mkdir -p $datadir/check
    bwa index $projdir/refs/ecoli_withplas.fasta.gz
    
    for i in NT278 NT279 NT296 NT297 NT298 NT299 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/ecoli_withplas.fasta.gz \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_bwa_ecoliplas.sorted.bam
	samtools index $datadir/check/${i}_bwa_ecoliplas.sorted.bam
    done
fi

if [ $1 == bowtie_withplas ] ; then
    bowtie2-build $projdir/refs/ecoli_withplas.fasta.gz $projdir/refs/ecoli_withplas
    for i in NT278 NT279 NT296 NT297 NT298 NT299 ;
    do
	bowtie2 \
	    -p 36 \
	    -x $projdir/refs/ecoli_withplas \
	    -1 $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    -2 $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_bt2_ecoliplas.sorted.bam
	samtools index $datadir/check/${i}_bt2.sorted_ecoliplas.bam
    done
fi	


if [ $1 == align_plas ] ; then
    mkdir -p $datadir/check
    bwa index $projdir/refs/ptkei-dest-cre-loxp-sfgfp.fasta.gz
    
    for i in NT278 NT279 NT296 NT297 NT298 NT299 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/ptkei-dest-cre-loxp-sfgfp.fasta.gz \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_bwa_plas.sorted.bam
	samtools index $datadir/check/${i}_bwa_plas.sorted.bam
    done
fi


if [ $1 == construct_plas ] ; then
    bwa index $projdir/refs/construct1_plas.fa.gz
    for i in NT278 NT296 NT297 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/construct1_plas.fa.gz \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_constructplas.sorted.bam
	samtools index $datadir/check/${i}_constructplas.sorted.bam
    done

    bwa index $projdir/refs/construct2_plas.fa.gz
    for i in NT279 NT298 NT299 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/construct2_plas.fa.gz \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_constructplas.sorted.bam
	samtools index $datadir/check/${i}_constructplas.sorted.bam
    done
fi

if [ $1 == ins_plas ] ; then
    bwa index $projdir/refs/ins1_plas.fa
    for i in NT278 NT296 NT297 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/ins1_plas.fa \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_insplas.sorted.bam
	samtools index $datadir/check/${i}_insplas.sorted.bam
    done

    bwa index $projdir/refs/ins2_plas.fa
    for i in NT279 NT298 NT299 ;
    do
	bwa mem \
	    -t 36 \
	    $projdir/refs/ins2_plas.fa \
	    $datadir/trimmed/${i}_fwd_paired.fq.gz \
	    $datadir/trimmed/${i}_rev_paired.fq.gz | \
	    samtools view -@ 36 -b | \
	    samtools sort -@ 36 -o $datadir/check/${i}_insplas.sorted.bam
	samtools index $datadir/check/${i}_insplas.sorted.bam
    done
fi
