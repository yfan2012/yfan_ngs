library(tidyverse)
library(RColorBrewer)


datadir='/dilithium/Data/NGS/projects/dunlop_rna/petriseq'
dbxdir='~/Dropbox/yfan/dunlop/single_cell/petriseq'

gene_entropy <- function(geneinfo) {
    interval=5
    high=max(geneinfo)
    expranges=tibble(low=seq(0, high, interval),
                     high=seq(interval, high+interval, interval))
    exp=expranges %>%
        rowwise() %>%
        mutate(count=sum(geneinfo>=low & geneinfo<high)) %>%
        ungroup() %>%
        mutate(frac=count/sum(count))
    exp$frac[exp$frac==0]=1 ##since p*log(p)=0 where p=0 by convention
    exp=exp %>%
        mutate(logfrac=log(frac, 10)) %>%
        mutate(ent=-frac*logfrac)

    entropy=sum(exp$ent)
    return(entropy)
}

cell_entropy <- function(nk) {
    tot=rowSums(nk %>% select(-cell))
    tot[tot==0]=1 #to avoid nans later
    
    ##calculate pi
    nk=nk %>% column_to_rownames('cell')
    p=nk/tot

    ##since p*log(p)=0 where p=0, by convention
    p[p==0]=1
    
    ##elementwise multiplication
    ent=as_tibble(enframe(rowSums(-p*log(p))))
    return(ent)
}

per_gene_info <- function(nk) {
    ##get a bunch of info per gene
    ##need nk matrix without cells column
    avg=nk %>%
        summarise(across(everything(), mean)) %>%
        gather()
    names(avg)=c('name', 'avg')
    std=nk %>%
        summarise(across(everything(), var)) %>%
        gather()
    names(std)=c('name', 'variance')
    ent=nk %>%
        summarise(across(everything(), gene_entropy)) %>%
        gather()
    names(ent)=c('name', 'ent')
    
    genes=as_tibble(enframe(colSums(nk!=0))) %>%
        rename(numcells=value) %>%
        mutate(numcells=numcells+1) %>%
        full_join(avg, by='name') %>%
        full_join(std, by='name') %>%
        full_join(ent, by='name') %>%
        mutate(rv=variance/avg)
    return(genes)
}    

per_cell_info <- function(nk) {
    ##get a bunch of info per cell
    ##need nk matrix with cells column
    ent=cell_entropy(nk) %>%
        rename(ent=value) %>%
        rename(cell=name)
    cells=nk %>%
        gather('gene', 'counts', -cell) %>%
        group_by(cell) %>%
        summarise(numgenes=sum(counts!=0),
                  numcounts=sum(counts), 
                  avg=mean(counts),
                  variance=var(counts)) %>%
        mutate(rv=variance/avg) %>%
        full_join(ent, by='cell')
    return(cells)
}

nonrrna <- function(nk) {
    ##get total counts per cell that are non-rrna
    notribo=enframe(rowSums(nk %>%
                            select(-c('cell', 'U00096:rRNA', 'CP000255:rRNA'))))
    return(notribo)
}
    
    
##check out supp data nk matrix, and see how many times 
exps=c('species_mix', 'growth_mix', 'growth_light_mix')
genes_of_interest=c('gadAXW', 'araC', 'recAX', 'rpoH')

cells_per_gene=tibble(name=as.character(),
                      numcells=as.integer(),
                      avg=as.numeric(),
                      variance=as.numeric(),
                      ent=as.numeric(),
                      rv=as.numeric(),
                      experiment=as.character())

genes_per_cell=tibble(cell=as.character(),
                      numgenes=as.integer(),
                      numcounts=as.integer(),
                      avg=as.numeric(),
                      variance=as.numeric(),
                      ent=as.numeric(),
                      rv=as.numeric(),
                      experiment=as.character())
geneinfo=tibble(name=as.character(),
              count=as.integer(),
              experiment=as.character())
allcounts=tibble(name=as.character(),
              count=as.integer(),
              experiment=as.character())
allnonribo=tibble(name=as.integer(),
                  value=as.numeric(), 
                  experiment=as.character())
for (i in exps) {
    nkfile=file.path(datadir, 'supp_data', paste0(i, '.tsv'))
    nk=read_table2(nkfile)
    nknames=c('cell', names(nk))
    names(nk)=nknames
    
    genes=per_gene_info(nk %>% select(-cell)) %>%
        mutate(experiment=i) 
    cells_per_gene=bind_rows(cells_per_gene, genes)

    genescsv=file.path(dbxdir, paste0(i,'.genes.csv'))
    write_csv(genes, genescsv)
    
    cells=per_cell_info(nk) %>%
        mutate(experiment=i)
    genes_per_cell=bind_rows(genes_per_cell, cells)

    interest=which(rowSums(outer(nknames, genes_of_interest, str_detect))>0)
    genedf=gather(nk[,interest]) %>%
        rename('name'='key', 'count'='value') %>%
        mutate(experiment=i)
    geneinfo=bind_rows(geneinfo, genedf)

    expcounts=gather(nk %>% select(-cell)) %>%
        rename('name'='key', 'count'='value') %>%
        mutate(experiment=i)
    allcounts=bind_rows(allcounts, expcounts)

    notribo=nonrrna(nk) %>%
        mutate(experiment=i)
    allnonribo=bind_rows(allnonribo, notribo)
}


pergenefile=file.path(dbxdir, 'cells_per_gene.pdf')
pdf(pergenefile, h=7, w=12)
ggplot(cells_per_gene, aes(x=numcells, colour=experiment, fill=experiment, alpha=.25)) +
    geom_density() +
    ggtitle('Cells per gene') +
    xlab('Number of Cells') +
    scale_x_log10() +
    scale_fill_brewer(palette = 'Set2') +
    scale_colour_brewer(palette = 'Set2') +
    theme_bw()
dev.off()

geneinfo=geneinfo %>%
    mutate(pseudocount=count+1)
goifile=file.path(dbxdir, 'genes_of_interest.pdf')
pdf(goifile, h=7, w=12)
for (i in exps) {
    plot=ggplot(geneinfo %>% filter(experiment==i), aes(x=pseudocount, colour=name, fill=name, alpha=.15)) +
        geom_density() +
        ggtitle(paste0('Expermiment:', i)) +
        xlab('Expression Counts') +
        scale_x_log10() +
        scale_fill_brewer(palette = 'Set2') +
        scale_colour_brewer(palette = 'Set2') +
        theme_bw()
    print(plot)
}
dev.off()



goitabfile=file.path(dbxdir, 'genes_of_interest.csv')
len=length(table(geneinfo$count))
histinfo=tibble(expression_count=as.character(),
                numcells=as.integer(),
                experiment=as.character(),
                gene=as.character())

for (i in exps) {
    genes=unique(geneinfo$name)
    for (j in genes){
        counthist=enframe(table(geneinfo %>% filter(experiment==i, name==j) %>% select(count))) %>%
            mutate(experiment=i) %>%
            mutate(gene=j) %>%
            mutate(value=as.integer(value))
        names(counthist)=c('expression_count','numcells', 'experiment', 'gene')
        histinfo=bind_rows(histinfo, counthist)
    }
}
histord=histinfo[, c(3,4,1,2)]
write_csv(histord, goitabfile)




highestcsv=file.path(dbxdir, 'highest_expression.csv')
highest=allcounts %>%
    group_by(experiment) %>%
    arrange(-count) %>%
    slice(1:5)
write_csv(highest, highestcsv)
highest_genescsv=file.path(dbxdir, 'highest_genes.csv')
highest_genes=allcounts %>%
    group_by(experiment, name) %>%
    summarise(maxcount=max(count)) %>%
    arrange(-maxcount) %>%
    slice(1:5)
write_csv(highest_genes, highest_genescsv)

allnonribo=allnonribo %>%
    mutate(pseudocount=value+1)
percellcountspdf=file.path(dbxdir, 'counts_per_cell.pdf')
pdf(percellcountspdf, w=12, h=7)
ggplot(allnonribo, aes(x=pseudocount, colour=experiment, fill=experiment, alpha=.25)) +
    geom_density() +
    ggtitle('Non-ribo counts per cell') +
    xlab('Pseudocounts') +
    scale_x_log10() +
    scale_fill_brewer(palette = 'Set2') +
    scale_colour_brewer(palette = 'Set2') +
    theme_bw()
dev.off()

    
genepos=which(rowSums(outer(cells_per_gene$name, genes_of_interest, str_detect))>0)
goi=cells_per_gene[genepos,] %>%
    rowwise() %>%
    mutate(rank=sum(avg<=cells_per_gene$avg[cells_per_gene$experiment==experiment]))
goirankcsv=file.path(dbxdir, 'sc_rank.csv')
write_csv(goi, goirankcsv)



###Get info on specific genes
genes=c('U00096:fliAZ-tcyJ','U00096:fliFGHIJK', 'U00096:fliLMNOPQR', 'U00096:fliC', 'U00096:motAB-cheAW', 'U00096:tar-tap-cheRBYZ')
i=exps[2]
nkfile=file.path(datadir, 'supp_data', paste0(i, '.tsv'))
nk=read_table2(nkfile)
nknames=c('cell', names(nk))
names(nk)=nknames

genecols=which(nknames %in% genes)
geneinfo=nk[,genecols]
growthgoicsv=file.path(dbxdir, 'growth_mix_goi.csv')
write_csv(geneinfo, growthgoicsv)

expcells=geneinfo[which(rowSums(geneinfo)!=0),]
patcells=expcells[which(rowSums(expcells)>1),]
                  
    

###get info on control genes
genes=c('U00096:elaB', 'U00096:phoPQ', 'U00096:lldPRD', 'U00096:wza-wzb-wzc-wcaAB', 'U00096:csgDEFG', 'U00096:norVW')
genecols=which(nknames %in% genes)
geneinfo=nk[,genecols]
growthgoicsv=file.path(dbxdir, 'growth_mix_goi_control.csv')
write_csv(geneinfo, growthgoicsv)

