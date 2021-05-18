library(tidyverse)
library(ShortRead)
library(cowplot)

dbxdir='~/Dropbox/yfan/dunlop/insert'
rawdir='/mithril/Data/NGS/Raw/210510_dunlop_nathan'
datadir='/mithril/Data/NGS/projects/dunlop_insert/run2/positions'
samps=c('NT284', 'NT285', 'NT286', 'NT287')

poscols=c('readname', 'pair', 'crepos', 'orient', 'redund', 'insert')


##filter for reads that actually have the interface
for (i in samps) {
    posfile=file.path(datadir, paste0(i, '.positions.csv'))
    posdata=read_csv(posfile, col_names=poscols) %>%
        filter(crepos!='None') %>%
        filter(redund=='FALSE')
    filtposfile=file.path(datadir, paste0(i, '.crepos_filt.positions.csv'))
    write_csv(posdata, filtposfile)
}
for (i in samps) {
    posfile=file.path(datadir, paste0(i, '.positions.csv'))
    posdata=read_csv(posfile, col_names=poscols) %>%
        filter(crepos!='None') %>%
        filter(orient=='True') %>%
        filter(redund=='FALSE')
    filtposfile=file.path(datadir, paste0(i, '.filtered.positions.csv'))
    write_csv(posdata, filtposfile)
}


####look at coverage
datadir='/mithril/Data/NGS/projects/dunlop_insert/run2/cov'
cov_cols=c('chrname', 'pos', 'cov')
covplots=list()
for (i in samps) {
    insertcovfile=file.path(datadir, paste0(i, '.genomecov'))
    ##plascovfile=file.path(datadir, paste0(i, 'plasmid.genomecov'))
    insert=read_tsv(insertcovfile, col_names=cov_cols) %>%
        mutate(samp=i)
    plot=ggplot(insert, aes(x=cov, colour=chrname, fill=chrname, alpha=.3)) +
        geom_density() +
        ggtitle(i) +
        scale_colour_brewer(palette='Set2') +
        scale_fill_brewer(palette='Set2') +
        theme_bw()
    covplots[[i]]=plot
}
covplotspdf=file.path(dbxdir, 'covplots.pdf')
pdf(covplotspdf, h=5, w=11)
plot_grid(covplots[[1]], covplots[[2]], covplots[[3]], covplots[[4]], ncol=2)
dev.off()


####get crepos 41 and 46 as an example
samp='NT284'
posfile=file.path(datadir, 'NT284.positions.csv')
posdata1=read_csv(posfile, col_names=poscols) %>%
    filter(crepos==41 | crepos==46) %>%
    filter(pair=='read1') %>%
    filter(redund=='FALSE')
posdata2=read_csv(posfile, col_names=poscols) %>%
    filter(crepos==41 | crepos==46) %>%
    filter(pair=='read2') %>%
    filter(redund=='FALSE')

fqfile1=file.path(rawdir, 'NT284_S1_L001_R1_001.fastq.gz')
fqfile2=file.path(rawdir, 'NT284_S1_L001_R2_001.fastq.gz')
fq1=readFastq(fqfile1)
fq2=readFastq(fqfile2)

read1=sapply(strsplit(as.character(id(fq1)), split=' ', fixed=TRUE), head, 1)
read2=sapply(strsplit(as.character(id(fq2)), split=' ', fixed=TRUE), head, 1)

wantr1=which(read1 %in% posdata1$readname)
wantr2=which(read2 %in% posdata2$readname)

names=c(id(fq1[wantr1]), id(fq2[wantr2]))
reads=c(sread(fq1[wantr1]), sread(fq2[wantr2]))
names(reads)=names
writefile=file.path(datadir, 'NT284_example.41.46.fasta')
writeXStringSet(reads, writefile, format='fasta')



