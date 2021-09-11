library(tidyverse)

datadir='/mithril/Data/NGS/projects/dunlop_insert/run3/exact'
dbxdir='~/gdrive/dunlop/insert'

naivenames=c('NT278', 'NT279')
poscols=c('name', 'pos', 'side', 'strand', 'read')
naive=tibble()

for (i in naivenames) {
    posfile=file.path(datadir, paste0(i, '.positions.csv'))
    naiveinfo=read_csv(posfile, col_names=poscols) %>%
        mutate(samp=i) %>%
        mutate(aapos=((pos-2)/3)+1) %>%
        mutate(frame=aapos%%1==0)
    naive=bind_rows(naive, naiveinfo)
}

covinfo=naive %>%
    group_by(pos, samp) %>%
    summarise(numreads=n())


##plot coverage histogram
covhistpdf=file.path(dbxdir, 'covhist_naive.pdf')
pdf(covhistpdf, w=15, h=7)
ggplot(covinfo, aes(x=numreads, colour=samp, fill=samp, alpha=.2)) +
    geom_histogram() +
    facet_wrap(~samp) +
    theme_bw()
dev.off()


##amino acid space
aacovinfo=naive %>%
    filter(frame==TRUE) %>%
    group_by(aapos, samp) %>%
    summarise(numreads=n())

for (i in seq(10,322,1)) {
    posmag=aacovinfo %>%
        filter(aapos==i) %>%
        filter(samp=='NT278')
    if (dim(posmag)[1]==0) {
        aacovinfo=bind_rows(aacovinfo, tibble(aapos=i, samp='NT278', numreads=0))
    }
    posmag=aacovinfo %>%
        filter(aapos==i) %>%
        filter(samp=='NT279')
    if (dim(posmag)[1]==0) {
        aacovinfo=bind_rows(aacovinfo, tibble(aapos=i, samp='NT279', numreads=0))
    }
}

covhistaapdf=file.path(dbxdir, 'covhist_naive_aa.pdf')
pdf(covhistaapdf, w=15, h=7)
ggplot(aacovinfo, aes(x=numreads, colour=samp, fill=samp, alpha=.2)) +
    geom_histogram() +
    facet_wrap(~samp) +
    theme_bw()
dev.off()
