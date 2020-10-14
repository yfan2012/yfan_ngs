library(tidyverse)

datadir='/mithril/Data/NGS/projects/dunlop_insert/positions/'
outdir='~/Dropbox/yfan/dunlop/insert/'
samps=c('B1_S1', 'B2_S2', 'B3_S3', 'B4_S4')


cols=c('readname', 'pair', 'crepos', 'orientmatch', 'pair_redundant')
allpos=tibble(readname=as.character(),
              pair=as.character(),
              crepos=as.numeric(),
              orientmatch=as.logical(),
              pair_redundant=as.logical(),
              samp=as.character())

for (i in samps) {
    posfile=paste0(datadir, i, '.positions.csv')
    pos=read_csv(posfile, col_names=cols) %>%
        filter(crepos!='None') %>%
        mutate(crepos=as.numeric(crepos)) %>%
        mutate(orientmatch=as.logical(orientmatch)) %>%
        mutate(samp=i)

    allpos=rbind(allpos, pos)
}

outfile=paste0(outdir, 'cre_insertion_locs.pdf')
pdf(outfile, h=7, w=13)
plot=ggplot(allpos, aes(x=crepos, colour=samp, fill=samp, alpha=.2)) +
    geom_density() +
    xlab('cre position') +
    ggtitle('cre insert locations') +
    theme_bw()
print(plot)
dev.off()
