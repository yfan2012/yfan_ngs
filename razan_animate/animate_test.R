library(ggplot2)
library(tidyverse)
library(gganimate)
theme_set(theme_bw())

dbxdir='~/gdrive/dunlop/razan_animation'
deadfile=file.path(dbxdir, 'gadXrecA-combdeadcells.csv')
livefile=file.path(dbxdir, 'gadXrecA-comblivecells.csv')

deaddata=read_csv(deadfile, col_names=F)
livedata=read_csv(livefile, col_names=F)
##separation at column SS and AMK, which is X513 and X1026

cfp=deaddata[,1:512]
yfp=deaddata[,513:1024]
len=deaddata[,1025:1536]

per.cell.data <- function(color1, color2, len, samp) {
    ##puts together long data - time series per cell
    percell=tibble(c1=as.numeric(),
                   c2=as.numeric(),
                   l=as.numeric(),
                   c=as.numeric())
    for (i in 1:dim(cfp)[1]) {
        cell=tibble(c1=as.numeric(color1[i,]), 
                    c2=as.numeric(color2[i,]),
                    l=as.numeric(len[i,]),
                    c=i)
        percell=bind_rows(percell, cell)
    }

    percell=percell %>%
        mutate(rec=case_when(is.na(l) ~ .2,
                             TRUE ~ .8)) %>%
        mutate(samp=samp)

    ##deal with NA
    for (i in 1:dim(percell)[1]) {
        if (is.na(percell$l[i])) {
            percell$l[i]=percell$l[i-1]
            percell$c1[i]=percell$c1[i-1]
            percell$c2[i]=percell$c2[i-1]
        }
    }
    return(percell)
}


percell=per.cell.data(livedata[,1:512], livedata[,513:1024], livedata[1025:1536], 'live')
percell=per.cell.data(deaddata[,1:512], deaddata[,513:1024], deaddata[1025:1536], 'dead')
