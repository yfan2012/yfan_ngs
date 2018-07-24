library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)

vcf='/home/yfan/Dropbox/yfan/sears_fmt/parsnp_ecoli_new/ecoli.vcf'

snptab=read_tsv(vcf, comment='##') %>%
    filter(FILTER=='PASS')


snps=matrix(nrow=17, ncol=17)
for (i in 1:17) {
    for (j in 1:17) {
        snps[i, j]=sum(snptab[,i+9]!=snptab[,j+9])
    }
}

names=colnames(snptab)[10:26]

for (i in 1:length(names)){
    newname=gsub('.spades.fasta', '',names[i])
    names[i]=newname
}
names[1]='NC101'

rownames(snps)=names
colnames(snps)=names

pdf('/home/yfan/Dropbox/yfan/sears_fmt/numsnps.pdf', width=20, height=7)
grid.table(snps)
dev.off()

