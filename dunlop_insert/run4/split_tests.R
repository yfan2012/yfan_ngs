library(tidyverse)

dbxdir='~/gdrive/dunlop/insert/plots'
datadir='/mithril/Data/NGS/projects/dunlop_insert/run4/exact'
samps=c('NT350.positions.csv','NT351.positions.csv','NT358.positions.csv','NT359.positions.csv')

exactcols=c('readname', 'split', 'side', 'orient', 'read')

splits=NULL
for (i in samps) {
    name=strsplit(i, '.', fixed=TRUE)[[1]][1]
    sampfile=file.path(datadir, i)
    splitinfo=read_csv(sampfile, col_names=exactcols) %>%
        mutate(samp=name) %>%
        select(-readname, -side, -orient, -read)
    splits=bind_rows(splits, splitinfo)
}

splitcounts=splits %>%
    group_by(samp, split) %>%
    summarise(count=n())

sums=splitcounts %>%
    group_by(samp) %>%
    summarise(allcounts=sum(count))
adjustscale=sums$allcounts[3]/sums$allcounts[4]
lightcounts=splitcounts %>%
    filter(samp=='NT358') %>%
    rename(light=count) %>%
    ungroup() %>%
    select(-samp)
darkcounts=splitcounts %>%

    filter(samp=='NT359') %>%
    mutate(dark=count*adjustscale) %>%
    select(-count) %>%
    ungroup() %>%
    select(-samp)

allcounts=full_join(darkcounts, lightcounts, by=c('split'))
allcounts[is.na(allcounts)]=0

countsfrac=allcounts %>%
    mutate(total=light+dark) %>%
    mutate(fraclight=light/total) %>%
    arrange(-fraclight)

countspdf=file.path(dbxdir, 'light_fraction.pdf')
pdf(countspdf, h=3, w=50)
plot=ggplot(countsfrac, aes(x=split, y=fraclight, alpha=.2)) +
    geom_bar(stat='identity') +
    xlab('split (nt)') +
    ylab('light fraction') +
    theme_bw()
print(plot)
dev.off()

##recall that 3*aa-1=nt
