library(Matrix)
library(tidyverse)

dbxdir='~/gdrive/cho'
rawdir='/atium/Data/projects/ambic_cho/220519_cho_igg/g1ks'
mitodir='/atium/Data/projects/ambic_cho/220503_alice_cho_sc/cho_k1gshd_analysis/'
mitofile=file.path(mitodir,'220519_cho_day90','mito_genes.txt')
matrix0dir=file.path(rawdir, 'CHO_day0/outs/filtered_feature_bc_matrix')
matrix9dir=file.path(rawdir, 'CHO_day90/outs/filtered_feature_bc_matrix')


testdir='/atium/Data/projects/ambic_cho_yfan/cellranger'
##matrix0dir=file.path(testdir, 'CHO_day0_count/outs/filtered_feature_bc_matrix')
##matrix90dir=file.path(testdir, 'CHO_day90_count/outs/filtered_feature_bc_matrix')

##check out the matrix first
##stolen from https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices#r-load-mat

mitodata=read_csv(mitofile, col_names=F)
colnames(mitodata)='mito'
mitogenes=as.vector(mitodata$mito)

loadmats <- function(matrixdir) {
    barcode.path=paste0(matrixdir, "/barcodes.tsv.gz")
    features.path=paste0(matrixdir, "/features.tsv.gz")
    matrix.path=paste0(matrixdir, "/matrix.mtx.gz")
    
    mat=readMM(file=matrix.path)
    feature.names = read.delim(features.path,
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path,
                               header = FALSE,
                               stringsAsFactors = FALSE)
    colnames(mat)=barcode.names$V1
    rownames(mat)=feature.names$V2
    
    return(mat)
}

mat0=loadmats(matrix0dir)
mat9=loadmats(matrix9dir)


genes.per.cell <- function(mat) {
    posgenes=mat!=0
    gpc=colSums(posgenes)
    return(gpc)
}

tscripts.per.cell <- function(mat) {
    tpc=colSums(mat)
    return(tpc)
}


mito.per.cell <- function(mat, mitogenes) {
    ##percent mitochondria counts per cell
    tpc=colSums(mat)
    mitomat=mat[rownames(mat) %in% mitogenes,]
    mpc=colSums(mitomat)

    pmpc=mpc/tpc
    return(pmpc)
}

filter.cells <- function(mat, mitogenes, gene.cut, mito.cut) {
    ##filter cells - gene.cut is c(min, max) for num genes, mito.cut is frac of mito
    ##mito.cut is a fraction - .05 for 5% cutoff
    gpc=genes.per.cell(mat)
    keepgpc=gpc[gpc>gene.cut[1] & gpc<gene.cut[2]]

    pmpc=mito.per.cell(mat, mitogenes)
    keeppmpc=pmpc[pmpc<mito.cut]

    keepcells=keeppmpc[names(keeppmpc) %in% names(keepgpc)]

    filtmat=mat[, keepcells]
    return(filtmat)
}


gene.cut=c(1500, 7000)
mito.cut=.1
filt0=filter.cells(mat0, mitogenes, gene.cut, mito.cut)
filt9=filter.cells(mat9, mitogenes, gene.cut, mito.cut)



####plotting
percell0=tibble(gpc=genes.per.cell(mat0), 
                tpc=tscripts.per.cell(mat0), 
                pmpc=mito.per.cell(mat0, mitogenes),
                label='day0')
percell9=tibble(gpc=genes.per.cell(mat90), 
                tpc=tscripts.per.cell(mat90), 
                pmpc=mito.per.cell(mat90, mitogenes),
                label='day90')
percell=bind_rows(percell0, percell9)


controlplotspdf=file.path(dbxdir, 'controlplots.pdf')
pdf(controlplotspdf, h=8, w=11)
gpcplot=ggplot(percell, aes(x=label, y=gpc, colour=label, fill=label, alpha=.2)) +
    geom_violin(alpha = 0.5) +
    geom_jitter(position = position_jitter(seed = 1, width = 0.2), size=.3) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    ggtitle('genes per cell') +
    theme_bw()
plot(gpcplot)
tpcplot=ggplot(percell, aes(x=label, y=tpc, colour=label, fill=label, alpha=.7)) +
    geom_violin(alpha = 0.5) +
    geom_jitter(position = position_jitter(seed = 1, width = 0.2), size=.3) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    ggtitle('transcripts per cell') +
    theme_bw()
plot(tpcplot)
mtplot=ggplot(percell, aes(x=label, y=pmpc, colour=label, fill=label, alpha=.7)) +
    geom_violin(alpha = 0.5) +
    geom_jitter(position = position_jitter(seed = 1, width = 0.2), size=.3) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    ggtitle('percent mito per cell') +
    theme_bw()
plot(mtplot)
dev.off()







