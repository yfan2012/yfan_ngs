library(ggplot2)
library(tidyverse)
library(gganimate)
library(gifski)
##library(av)
theme_set(theme_bw())

dbxdir='~/gdrive/dunlop/razan_animation'
deadfile=file.path(dbxdir, 'gadXrecA-combdeadcells.csv')
livefile=file.path(dbxdir, 'gadXrecA-comblivecells.csv')

deaddata=read_csv(deadfile, col_names=F)
livedata=read_csv(livefile, col_names=F)
##separation at column SS and AMK, which is X513 and X1026

per.cell.data <- function(color1, color2, len, samp) {
    ##puts together long data - time series per cell
    percell=tibble(c1=as.numeric(),
                   c2=as.numeric(),
                   l=as.numeric(),
                   c=as.numeric(),
                   time=as.numeric())
    for (i in 1:dim(color1)[1]) {
        cell=tibble(c1=as.numeric(color1[i,]), 
                    c2=as.numeric(color2[i,]),
                    l=as.numeric(len[i,]),
                    c=i)
        cell$time=c(1:dim(cell)[1])
        percell=bind_rows(percell, cell)
    }

    percell=percell %>%
        mutate(rec=case_when(is.na(l) ~ .15,
                             TRUE ~ .6)) %>%
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


livepercell=per.cell.data(livedata[,1:512], livedata[,513:1024], livedata[1025:1536], 'live')
deadpercell=per.cell.data(deaddata[,1:512], deaddata[,513:1024], deaddata[1025:1536], 'dead')

allcells=bind_rows(livepercell, deadpercell)

plot=ggplot(allcells, aes(x=c1, y=c2, size=l, colour=samp, alpha=rec))+
    geom_point() +
    scale_colour_brewer(palette='Set2') +
    scale_y_log10() +
    scale_x_log10() +
    scale_alpha_continuous(range = c(.15, 0.6)) +
    labs(x = "Color1", y = "Color2") +
    transition_time(time) +
    labs(title = "Frame: {frame_time}")

animpath=file.path(dbxdir, 'test.gif')
future::plan("multiprocess", workers = 36L)
anim_save(animpath, plot, nframes=500, fps=10, renderer=gifski_renderer())

                                              
