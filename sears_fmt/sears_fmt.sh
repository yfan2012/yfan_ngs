#!/bin/bash


##Set up file locations
rawdir=/dilithium/Data/NGS/Raw/171115_sears_fmt
datadir=/dilithium/Data/NGS/Aligned/171115_sears_fmt
refpre=/mithril/Data/NGS/Reference/ecoli/ecoli
ref=/mithril/Data/NGS/Reference/ecoli/ecoli.fasta

##Run the pipeline for each sample. Actually run the loop over only the read1 files to avoid running twice for each sample. 
for i in $rawdir/*R1_001.fastq.gz;
do
    ##right now, the variable i is a string that looks like this: '/dilithium/Data/NGS/Raw/171115_sears_fmt/1D12_S1_L001_R1_001.fastq.gz'
    ##I just want the part that says '1D12_S1_L001' to use as a sample label for the output files that come out later, so I'm making a new variable called prefix, which contains that label. 
    ## echo ${i#$rawdir/} ---> cuts off everything in rawdir(above), so the label looks like this: 1D12_S1_L001_R1_001.fastq.gz
    ## cut -d '_' -f 1,2 ---> cuts '1D12_S1_L001_R1_001.fastq.gz' wherever there's an underscore, and keeps the first and second fields, so now the label looks like this: 1D12_S1

    prefix=`echo ${i#$rawdir/} | cut -d '_' -f 1,2`
    echo the full path of the file being worked on by the for loop: $i
    echo just the name of the file: ${i#$rawdir/}
    echo just the label we get from the file name: $prefix
    

    ##If the data directory doesn't exist, make it.
    mkdir -p $datadir
   

    
    ##Read mapping/alignment
    ##An 'if statement' only runs the code if what's 'in' the statement is true.
    # the ! generally means 'not'
    # -f means 'file'
    ##So this this statement says: if there is no file called $datadir/$prefix.sorted.bam, then run bowtie2
    ##this is helpful for debugging. If you find problems with lines of code after this one, you can re-run the whole script, and it will automatically skip this step if it worked the first time. 

    if [ ! -f $datadir/$prefix.sorted.bam ] ; then
	##Run bowtie2, which gives you a sam file. Pipe (the '|' symbol) that sam file to 'samtools view' to convert it to a bam (which is a conpressed sam), and pipe  that bam to 'samtools sort.' 
	bowtie2 -p 12 -x $refpre -1 $rawdir/${prefix}_L001_R1_001.fastq.gz -2 $rawdir/${prefix}_L001_R2_001.fastq.gz | samtools view -bS | samtools sort -o $datadir/$prefix.sorted.bam
	samtools index $datadir/$prefix.sorted.bam
    fi


    ##Call variants
    if [ ! -f $datadir/mpileup_grubii/$prefix.cons.fq ] ; then
	##samtools mpileup: look at alignment and figure out coverage, overlaps, etc
	##bcftools call: based on mpileup info, make a decision about which positions are different between the sample reads and the reference genome
	##vcfutils.pl vcf2fq: convert vcf output of bcftools into fastq format
	samtools mpileup -uf $ref $datadir/$prefix.sorted.bam | bcftools call -c - | vcfutils.pl vcf2fq > $datadir/$prefix.cons.fq

	##filter based on coverage depth of 100
	samtools mpileup -uf $ref $datadir/$prefix.sorted.bam | bcftools call -c - | vcfutils.pl varFilter -D100 > $datadir/$prefix.vcf

	##convert fastq to fasta
	seqtk seq -a $datadir/$prefix.cons.fq > $datadir/$prefix.fasta
    fi
    
done


