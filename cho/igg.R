library(Matrix)
library(tidyverse)
library(qlcMatrix)


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

    filtmat=mat[, names(keepcells)]
    return(filtmat)
}


gene.cut=c(1500, 7000)
mito.cut=.1
filt0=filter.cells(mat0, mitogenes, gene.cut, mito.cut)
filt9=filter.cells(mat9, mitogenes, gene.cut, mito.cut)

normalize.counts.old <- function(mat, scalefactor) {
    ##normalize counts across cells
    tpc=tscripts.per.cell(mat)
    ivec=c()
    j=c()
    x=c()
    ##takes ~20s per 100 cells - might come back and parallel this
    for (i in 0:length(tpc)-1) {
        cellindex=which(mat@j==i)
        genecounts=mat@x[cellindex]
        normgenecounts=genecounts*scalefactor/tpc[i+1]
        x=c(x, normgenecounts)
        j=c(j, mat@j[cellindex])
        ivec=c(ivec, mat@i[cellindex])
        if (i %% 100 == 0 ) {
            print(i)
        }
    }

    norm=sparseMatrix(i=ivec, j=j, x=x, dims=c(20825, 4950), repr='T')
    
}

normalize.counts <- function(mat, scalefactor) {
    ##normalize counts across cells
    matdf=tibble(i=mat@i, j=mat@j, x=mat@x) %>%
        group_by(j) %>%
        summarise(i=i, norm=x*scalefactor/sum(x))

    norm=sparseMatrix(i=matdf$i+1, j=matdf$j+1, x=matdf$norm, dims=dim(mat), repr='T')
    rownames(norm)=rownames(mat)
    colnames(norm)=colnames(mat)
    return(norm)
}


scalefactor=10000
norm0=normalize.counts(filt0, scalefactor)
norm9=normalize.counts(filt9, scalefactor)

goi=rbind(tibble(expr=norm0['igg_hc',], cond='0', gene='igg_hc'),
          tibble(expr=norm0['igg_lc',], cond='0', gene='igg_lc'),
          tibble(expr=norm9['igg_hc',], cond='9', gene='igg_hc'),
          tibble(expr=norm9['igg_lc',], cond='9', gene='igg_lc'),
          tibble(expr=norm0['Actb',], cond='0', gene='Actb'),
          tibble(expr=norm9['Actb',], cond='9', gene='Actb'),
          tibble(expr=norm0['Gss',], cond='0', gene='Gss'),
          tibble(expr=norm9['Gss',], cond='9', gene='Gss')) %>%
    group_by(gene) %>%
    summarise(normexpr=expr/max(expr), cond=cond)


gene.correlation.old <- function(gene, mat) {
    ##calculates correlation between gene and all other genes
    ##this is too slow
    library(foreach)
    library(plyr)
    library(doParallel)
    cl=makeCluster(5)
    registerDoParallel(cl, cores=5)
    exprgene=mat[gene,]
    
    corrdf=foreach(i=1:dim(mat)[1], .combine=rbind) %dopar% {
        print(i)
        g2=mat[i,]
        coeff=cor(exprgene, g2, method='spearman')
        genename=rownames(mat)[i]
        return(data.frame(gene=genename, corr=coeff))
    }

    stopCluster(cl)
    rm(cl)
    return(coeff)
}


##taken from https://saket-choudhary.me/blog/2022/03/10/fast-sparsespearman/
SparsifiedRanks2 <- function(X) {
    if (class(X)[1] != "dgCMatrix") {
        X <- as(object = X, Class = "dgCMatrix")
    }
    non_zeros_per_col <- diff(x = X@p)
    n_zeros_per_col <- nrow(x = X) - non_zeros_per_col
    offsets <- (n_zeros_per_col - 1) / 2
    x <- X@x
    ## split entries to columns
    col_lst <- split(x = x, f = rep.int(1:ncol(X), non_zeros_per_col))
    ## calculate sparsified ranks and do shifting
    sparsified_ranks <- unlist(x = lapply(X = seq_along(col_lst),
                                          FUN = function(i) rank(x = col_lst[[i]]) + offsets[i]))
    ## Create template rank matrix
    X.ranks <- X
    X.ranks@x <- sparsified_ranks
    return(X.ranks)
}

SparseSpearmanCor2 <- function(X, Y = NULL, cov = FALSE) {
    ## Get sparsified ranks
    rankX <- SparsifiedRanks2(X)
    if (is.null(Y)){
        ## Calculate pearson correlation on rank matrices
        return (corSparse(X=rankX, cov=cov))
    }
    rankY <- SparsifiedRanks2(Y)
    return(corSparse( X = rankX, Y = rankY, cov = cov))
}

corr0=SparseSpearmanCor2(t(norm0))
corr9=SparseSpearmanCor2(t(norm9))
rownames(corr0)=rownames(norm0)
colnames(corr0)=rownames(norm0)
rownames(corr9)=rownames(norm9)
colnames(corr9)=rownames(norm9)

grab.gene.corr <- function(name, corrmat, samp) {
    genecorr=tibble(corr=corrmat[,name], gene=rownames(corrmat)) %>%
        filter(!is.na(corr)) %>%
        arrange(desc(corr)) %>%
        mutate(order=rank(desc(corr))) %>%
        mutate(name=name) %>%
        mutate(samp=samp)
    return(genecorr)
}

plotcorrdata=rbind(grab.gene.corr('igg_hc', corr0, 'day0'),
                   grab.gene.corr('igg_hc', corr9, 'day9'),
                   grab.gene.corr('igg_lc', corr0, 'day0'),
                   grab.gene.corr('igg_lc', corr9, 'day9'),
                   grab.gene.corr('Gss', corr0, 'day0'),
                   grab.gene.corr('Gss', corr9, 'day9'),
                   grab.gene.corr('Actb', corr0, 'day0'),
                   grab.gene.corr('Actb', corr9, 'day9'))

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


iggplotspdf=file.path(dbxdir, 'goi_dists.pdf')
pdf(iggplotspdf, h=8, w=13)
iggplot=ggplot(goi, aes(x=gene, y=normexpr, colour=cond, fill=cond, alpha=.2)) +
    geom_violin(alpha = 0.5) +
    ##geom_jitter(position = position_jitter(seed = 1, width = 0.2), size=.3) +
    scale_colour_brewer(palette='Set2') +
    scale_fill_brewer(palette='Set2') +
    ggtitle('genes per cell') +
    theme_bw()
print(iggplot)
dev.off()


corrplotspdf=file.path(dbxdir, 'corr_plot.pdf')
pdf(corrplotspdf, h=8, w=13)
corrplot=ggplot(plotcorrdata, aes(x=order, y=corr)) +
    geom_point(alpha=.4, size=.6) +
    facet_grid(samp ~ name) +
    scale_colour_brewer(palette='Set2') +
    theme_bw()
plot(corrplot)
dev.off()

