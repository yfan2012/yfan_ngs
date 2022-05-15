library(tidyverse)

projdir='/mithril/Data/NGS/projects/dunlop_insert'
datadir=file.path(projdir, 'run4')
dbxdir='~/gdrive/dunlop/insert/plots'
ref1=file.path(projdir, 'refs/construct1.fa')

samps=c('NT350', 'NT351', 'NT358', 'NT359')
naive=c(samps[1], samps[2])
poscols=c('readname', 'position', 'side', 'strand', 'read')

plotfile=file.path(dbxdir, 'posfreq.pdf')
pdf(plotfile, h=7, w=20)
for (i in naive) {
    posfile=file.path(datadir, 'exact', paste0(i, '.positions.csv'))
    posinfo=read_csv(posfile, col_names=poscols)
    posfreq=posinfo %>%
        group_by(position) %>%
        summarise(freq=n())
    plot=ggplot(data=posfreq, aes(x=position, y=freq, alpha=.2)) +
        geom_bar(stat='identity') +
        ggtitle(i) +
        theme_bw()
    print(plot)
}
dev.off()
    
    
