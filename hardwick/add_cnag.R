library(tidyverse)

datadir='/kyber/Data/NGS/projects/190513_hardwick'
dbxdir='~/Dropbox/yfan/hardwick'

samps=c(197,178, 1694, 6341)
combs=t(combn(samps,2))


##prepgff
gfffile=file.path(datadir, 'reference/CRNE_H99.gff')
gffcols=c('chrom', 'db', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')
gff=read_tsv(gfffile, col_names=gffcols, comment='#')

##from https://stackoverflow.com/questions/43014782/how-to-get-the-nth-element-of-each-item-of-a-list-which-is-itself-a-vector-of-u
find_locus <- function(gene, gff){
    row=gff$attribute[grepl(gene, gff$attribute, fixed=TRUE)]

    if (length(row)[1]>=1) {
        loci=str_split(row, ';')
        tags=sapply(loci, `[`, 5)
        splittags=str_split(tags, '=')
        cnag=sapply(splittags, `[`, 2)
        singlecnags=unique(cnag)
        report=paste0(singlecnags, collapse='|')
        return(report)
    }else{
        return('none')
    }
}



for (i in 1:dim(combs)[1]) {
    comb=combs[i,]
    filename=paste0(as.character(comb[1]), 'vs', as.character(comb[2]), '.csv')
    filepath=file.path(dbxdir, filename)
    info=read_csv(filepath) %>%
        rowwise() %>%
        mutate(cnag=find_locus(genename, gff))
    newpath=file.path(dbxdir, paste0('cnag.', filename))
    write_csv(info, newpath)
}
