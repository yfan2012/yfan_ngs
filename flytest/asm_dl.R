library(tidyverse)

asmdir='/kyber/Data/NGS/projects/190725_flywtf/ref'
asmfile=paste0(asmdir, '/genomes_euks.csv')
asmtab=read_csv(asmfile)

for (i in asmtab$GenBank){
    filename=strsplit(i, split='/', fixed=TRUE)[[1]][10]
    system(paste0('wget ', i, '/', filename, '_genomic.fna.gz ', asmdir, '/', filename, '_genomic.fna.gz '))
}
