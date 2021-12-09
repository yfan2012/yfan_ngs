library(tidyverse)
library(cowplot)

projdir='/mithril/Data/NGS/projects/dunlop_insert'
dbxdir='~/gdrive/dunlop/insert/plots'

sampinfo=tibble(runs=c('run3', 'run3', 'run4', 'run4'),
                naivests=c('NT278', 'NT279', 'NT350', 'NT351'))
pos_cols=c('readname', 'pos', 'side', 'orient', 'pair')

positions=NULL
for (i in 1:dim(sampinfo)) {
    samp=sampinfo[i,]
    sampfile=file.path(projdir, samp$runs, 'exact', paste0(samp$naivests, '.positions.csv'))
    pos=read_csv(sampfile, col_names=pos_cols) %>%
        group_by(pos) %>%
        summarise(cov=n()) %>%
        mutate(strain=samp$naivests)
    positions=bind_rows(positions, pos)
}

combs=combn(sampinfo$naivests, 2)

corr_sample <- function(posdata){
    cors=NULL
    corrdata=posdata %>%
        mutate(area=x*y) %>%
        arrange(-area)
    for (i in 1:dim(corrdata)[1]) {
        plotdata=corrdata[i:dim(corrdata)[1],]
        p=cor(plotdata$x, plotdata$y)
        cors=bind_rows(cors, tibble(rank=i, cor=p))
    }
    return(cors)
}

corr_sample_sum <- function(posdata) {
    cors=NULL
    corrdata=posdata %>%
        mutate(area=x+y) %>%
        arrange(-area)
    for (i in 1:dim(corrdata)[1]) {
        plotdata=corrdata[i:dim(corrdata)[1],]
        p=cor(plotdata$x, plotdata$y)
        cors=bind_rows(cors, tibble(rank=i, cor=p))
    }
    return(cors)
}
    
plots=c()
rankplots=c()
corplots=c()
sumcorplots=c()
for (i in 1:dim(combs)[2]){
    xstrain=combs[1,i]
    ystrain=combs[2,i]

    xcov=positions %>%
        filter(strain==xstrain) %>%
        select(-strain) %>%
        rename(x=cov)
    
    ycov=positions %>%
        filter(strain==ystrain) %>%
        select(-strain) %>%
        rename(y=cov)

    posdata=left_join(xcov, ycov)
    posdata[is.na(posdata)]=0
    p=cor(posdata$x, posdata$y)

    plot=ggplot(posdata, aes(x=x, y=y, alpha=.05, size=.5)) +
        geom_point(alpha=.5, size=.1) +
        xlab(xstrain) +
        ylab(ystrain) +
        ggtitle(paste0('pearson=', as.character(p))) +
        theme_bw()
    plots[[i]]=plot
    
    rankdata=posdata %>%
        arrange(-x) %>%
        mutate(xrank=seq(1,n(),1)) %>%
        arrange(-y) %>%
        mutate(yrank=seq(1,n(),1))
    rankp=cor(rankdata$xrank, rankdata$yrank, method='spearman')
    
    rankplot=ggplot(rankdata, aes(x=xrank, y=yrank, alpha=.05, size=.5)) +
        geom_point(alpha=.5, size=.1) +
        xlab(xstrain) +
        ylab(ystrain) +
        ggtitle(paste0('spearman=', as.character(rankp))) +
        theme_bw()
    rankplots[[i]]=rankplot

    cors=corr_sample(posdata)
    corplot=ggplot(cors, aes(x=rank, y=cor)) +
        geom_point(alpha=.5, size=.1) +
        xlab('Number of top points missing') +
        ylab('Pearson coefficient') +
        ggtitle(paste0(xstrain, ' vs ', ystrain)) +
        theme_bw()
    corplots[[i]]=corplot

    sumcors=corr_sample_sum(posdata)
    sumcorplot=ggplot(sumcors, aes(x=rank, y=cor)) +
        geom_point(alpha=.5, size=.1) +
        xlab('Number of top points missing') +
        ylab('Pearson coefficient') +
        ggtitle(paste0(xstrain, ' vs ', ystrain)) +
        theme_bw()
    sumcorplots[[i]]=sumcorplot

}


correlationfile=file.path(dbxdir, 'naivelib_cov_correlations.pdf')
pdf(correlationfile, h=8, w=15)
print(plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], ncol=3, align='v'))
dev.off()

rankcorfile=file.path(dbxdir, 'naivelib_covrank_correlations.pdf')
pdf(rankcorfile, h=8, w=15)
print(plot_grid(rankplots[[1]], rankplots[[2]], rankplots[[3]], rankplots[[4]], rankplots[[5]], rankplots[[6]], ncol=3, align='v'))
dev.off()

corrplotfile=file.path(dbxdir, 'naivelib_corrplot.pdf')
pdf(corrplotfile, h=8, w=15)
print(plot_grid(corplots[[1]], corplots[[2]], corplots[[3]], corplots[[4]], corplots[[5]], corplots[[6]], ncol=3, align='v'))
dev.off()

corrplotfile=file.path(dbxdir, 'naivelib_sumcorrplot.pdf')
pdf(corrplotfile, h=8, w=15)
print(plot_grid(sumcorplots[[1]], sumcorplots[[2]], sumcorplots[[3]], sumcorplots[[4]], sumcorplots[[5]], sumcorplots[[6]], ncol=3, align='v'))
dev.off()
