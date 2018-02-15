##blast for pks island

pksref=/mithril/Data/NGS/Reference/ecoli/pks_island.fasta
datadir=/dilithium/Data/NGS/Aligned/171115_sears_fmt
readdir=/dilithium/Data/NGS/Raw/171115_sears_fmt


##make a blast database from the genes
if [ $1 == blast ] ; then

    mkdir -p $datadir/blastdb
    mkdir -p $datadir/pks_blast
    mkdir -p $datadir/reads
    
    for i in $readdir/*fastq.gz ;
    do
	samp=` basename $i _001.fastq.gz `
	seqtk seq -a $i > $datadir/reads/$samp.fasta
	makeblastdb -in $datadir/reads/$samp.fasta -out $datadir/blastdb/$samp.db -dbtype nucl
	blastn -query $pksref -db $datadir/blastdb/$samp.db -outfmt 7 -out $datadir/pks_blast/$samp.pks.tsv -num_threads 12
    done
fi


if [ $1 == assembly ] ; then

    for i in $datadir/assemblies/all/*.fasta ;
    do
	samp=`basename $i .fasta`
	makeblastdb -in $i -out $datadir/blastdb/$samp.db -dbtype nucl
	blastn -query $pksref -db $datadir/blastdb/$samp.db -outfmt 7 -out $datadir/pks_blast/$samp.pks.tsv -num_threads 12
    done
fi


