suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library('tidyverse'))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("plyr"))

snp_annotate <- function(vcffile, gfffile, prefix) {
    ##calculate how many snps are in annotated regions regions

    vcf=read_tsv(vcffile, comment='##') %>%
        filter(FILTER=='PASS')

    names=colnames(vcf)[10:dim(vcf)[2]]
    for (i in 1:length(names)){
        newname=gsub('.fasta', '',names[i])
        names[i]=newname
    }
    colnames(vcf)[10:dim(vcf)[2]]=names

    qcol=match(prefix, colnames(vcf))
    snplist=vcf[vcf[,qcol]==1,]

    gff=read_tsv(gfffile, comment='#', col_names=FALSE)
    colnames(gff)=c('contig', 'caller', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene')

    inregion=adply(snplist, 1, function(x) {
        tiggenes=gff[gff$contig==x[1],]
        if(dim(tiggenes)[1]>0){
            ingene=sum(x[2]>tiggenes$start && x[2]<tiggenes$end)
            x$ingene=ingene
        }else{
            x$ingene=0
        }
        return(x)
    })
    return(sum(inregion$ingene))
}


##coding region of 178, test for all annotations
vcf178='~/Dropbox/yfan/hardwick/parsnp_long/strain_long_snps.vcf'
gff178="/kyber/Data/NGS/projects/190513_hardwick/augustus/178_all.gff"
prefix="197_over20k"

snp_annotate(vcf178, gff178, prefix)
