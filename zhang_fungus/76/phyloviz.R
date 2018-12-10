library(tidyverse)
library(reshape2)
library(ape)
library(prodlim)

allelefile='~/Dropbox/yfan/fungus_zhang/fungus_76/mlst/phyloviz/mlst_alleles_species.csv'

calls=read_csv(allelefile) %>%
    mutate(allelenum=as.numeric(sapply(strsplit(allele, '_'), function(x) x[2]))) %>%
    select(-len) %>%
    select(-allele)
 
phyloviz=spread(calls, gene, allelenum, fill=0)

write_delim(phyloviz, '~/Dropbox/yfan/fungus_zhang/fungus_76/mlst/phyloviz/phyloviz.tsv', delim='\t')


distfile='~/Dropbox/yfan/fungus_zhang/fungus_76/mlst/phyloviz/distanceMatrix.tab'
dist=as.dist(read.table(distfile, sep='\t'))
clust=hclust(dist)

tree=as.phylo(clust)

write.tree(phy=tree, file='~/Dropbox/yfan/fungus_zhang/fungus_76/mlst/phyloviz/phyloviztree.newick')

##checking for identical mlst calls
for (i in 1:dim(phyloviz)[1]) {
    row=phyloviz[i,2:dim(phyloviz)[2]]
    match=row.match(row, phyloviz[,2:dim(phyloviz)[2]])
    if (match < i) {
        print(paste0(phyloviz$sample[i], ' ', phyloviz$sample[match]))
    }
}

