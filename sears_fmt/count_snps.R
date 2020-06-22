library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)

vcf='/home/yfan/Dropbox/yfan/sears_fmt/parsnp_ecoli_new/ecoli.vcf'

snptab=read_tsv(vcf, comment='##') %>%
    filter(FILTER=='PASS')


snps=matrix(nrow=15, ncol=15)
for (i in 1:15) {
    for (j in 1:15) {
        snps[i, j]=sum(snptab[,i+9]!=snptab[,j+9])
    }
}

names=colnames(snptab)[10:24]

for (i in 1:length(names)){
    ##newname=strsplit(gsub('.spades.fasta', '',names[i]), split= '_', fixed=TRUE)[[1]][1]
    newname=gsub('.spades.fasta', '',names[i])
    names[i]=newname
}
names[1]='NC101'


rownames(snps)=names
colnames(snps)=names

##remove the 'new' samples since we trust that they're right
snps=snps[c(1,2,3,4,6,7,8,9,13,14,15),c(1,2,3,4,6,7,8,9,13,14,15)]
names=names[c(1,2,3,4,6,7,8,9,13,14,15)]

for (i in 1:length(names)){
    newname=strsplit(names[i], split= '_', fixed=TRUE)[[1]][1]
    names[i]=newname
}

names[8]='2R12C'
rownames(snps)=names
colnames(snps)=names


##re-order
roword=snps[order(snps[,7]),]
snps=roword[,order(roword[3,])]
snps=snps[c(3,1,2,4,5,6,7,8,10,11,9), c(3,1,2,4,5,6,7,8,10,11,9)]


##rename from 12 to 9
ordered_names=rownames(snps)
for (i in 1:length(ordered_names)){
    newname1=gsub('12', '9', ordered_names[i])
    newname=gsub('R', 'P', newname1)
    ordered_names[i]=newname
}

rownames(snps)=ordered_names
colnames(snps)=ordered_names

snps[lower.tri(snps,diag=TRUE)]=''

write.table(snps, '/home/yfan/Dropbox/yfan/sears_fmt/numsnps.csv', sep=',', col.names=NA)

pdf('/home/yfan/Dropbox/yfan/sears_fmt/numsnps.pdf', width=20, height=7)
grid.table(snps)
dev.off()

