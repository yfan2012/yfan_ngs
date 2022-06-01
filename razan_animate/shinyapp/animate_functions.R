library(ggplot2)
library(tidyverse)
library(gganimate)
##library(gifski)
library(av)
theme_set(theme_bw())

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


sep.colors <- function(celldata, tf, status) {
    ##tf is total frames
    sepdata=per.cell.data(celldata[,1:tf], celldata[,tf+1:tf+tf], celldata[,tf+tf+1:tf+tf+tf], status)
    return(sepdata)
}
    
    
#livepercell=sep.colors(livedata[,1:512], livedata[,513:1024], livedata[1025:1536], 'live')
#deadpercell=per.cell.data(deaddata[,1:512], deaddata[,513:1024], deaddata[1025:1536], 'dead')


make_movie <- function(livepercell, deadpercell) {
    allcells=bind_rows(livepercell, deadpercell)
    
    plot=ggplot(allcells, aes(x=c1, y=c2, size=l, colour=samp, alpha=rec))+
        geom_point() +
        scale_colour_brewer(palette='Set2') +
        ##scale_y_log10() +
        ##scale_x_log10() +
        scale_alpha_continuous(range = c(.15, 0.6)) +
        theme(legend.position="none") +
        labs(x = "Color1", y = "Color2") +
        transition_time(time) +
        labs(title = "Frame: {frame_time}")
    
    library('gganimateparallel')
    future::plan("multiprocess", workers = 36L)
    
    a=animate(plot, duration=40,
              fps=50,
              height=5,
              width=6,
              units="in",
              res=200,
              renderer=av_renderer())
    return(a)
    ##anim_save(animpath, a)
}

