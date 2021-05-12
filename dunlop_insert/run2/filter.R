library(tidyverse)

datadir='/mithril/Data/NGS/projects/dunlop_insert/run2/positions'
samps=c('NT284', 'NT285', 'NT286', 'NT287')

poscols=c('readname', 'pair', 'crepos', 'orient', 'redund')

for (i in samps) {
    posfile=file.path(datadir, paste0(i, '.positions.csv'))
    posdata=read_csv(posfile, col_names=poscols) %>%
        filter(crepos!='None') %>%
        filter(orient=='True') %>%
        filter(redund=='FALSE') %>%
        mutate(samp=i)
    filtposfile=file.path(datadir, paste0(i, '.filtered.positions.csv'))
    write_csv(posdata, filtposfile)
}

    
