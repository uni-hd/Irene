dPCA <- function (meta, bed, data, minlen = 50, maxlen = 2e4, lambda = 0.2, control="Healthy", transform = TRUE, trunc=TRUE, 
    nPaired = 0, nTransform = 0, nColMeanCent = 0, nColStand = 0, nMColMeanCent = 1, 
    nMColStand = 0, dSNRCut = 5, nUsedPCAZ = 0, nUseRB = 0, dPeakFDRCut = 0.5) 
{
    dd <- get.norm.data(meta, bed, data, lambda=lambda, control=control, transform=transform, trunc=trunc)
    bed <- dd$bed
    groupId <- as.numeric(as.factor(meta$condition))
    datasetId <- as.numeric(as.factor(meta$factor))
    condId <- unique(groupId)
    nGroupNum <- length(condId)
    nDatasetNum <- length(unique(datasetId))
    nSampleNum <- nrow(meta)
    sampleName <- rep("", nSampleNum)
    repNum <- matrix(unlist(lapply(condId, function(i) {
        x <- datasetId[groupId == i]
        unlist(lapply(1:max(datasetId), function(j) {
            length(which((x == j)))
        }))
    })), nrow = nGroupNum, byrow = TRUE)
    nLociNum <- nrow(bed)
    d <- dPCA_main_impl(nGroupNum, nDatasetNum, nSampleNum, 
            nPaired, nLociNum, groupId, datasetId, 
            repNum, sampleName, bed[, 1], bed[, 2], bed[, 3], 
            dd$data, nTransform, nColMeanCent, nColStand, nMColMeanCent, 
            nMColStand, dSNRCut, nUsedPCAZ, nUseRB, dPeakFDRCut)
    d$proj <- matrix(unlist(lapply(1:4, function(i) d[[i]])), 
            nrow = 4, byrow = TRUE, dimnames = list(names(d)[1:4]))
    d$Dobs <- matrix(d$Dobs, ncol = nDatasetNum, byrow = TRUE)
    d$PC <- matrix(d$PC, ncol = nDatasetNum, byrow = TRUE)
    d$lambda <- dd$lambda
    d$bed <- bed
    d
}

get.norm.data <- function (meta, bed, data, minlen = 50, maxlen = 2e4, lambda = 0.2, control="Healthy", transform = FALSE, trunc=FALSE) {
    mk <- sort(unique(meta$factor))
    len <- abs(bed[, 3] - bed[, 2]) + 1
    ix <- len>minlen & len<maxlen
    bed <- bed[ix, ]
    data <- data[ix, ] * 1e3/len[ix]
    z <- min(data[!is.na(data)])
    if (z <= 0) 
        data <- data - z + 1e-5
    if (lambda=="auto")
        lambda = get.lambda(data)
    if (transform)
        data <- boxCox(data, lambda)
    if (trunc)
        data[is.na(data)] = 0
    data <- normalizeQuantiles(data)
    list(bed=bed, data=data, lambda=lambda, Dobs=do.call(cbind,lapply(mk, function(i) 
        rowMeans(data[,meta$condition!=control & meta$factor==i]) - rowMeans(data[,meta$condition==control & meta$factor==i]))))
}

boxCox <- function(x,lambda=0.2) {
    (x^lambda - 1)/lambda
}

invboxCox <- function(x,lambda=0.2) {
    exp(log(abs(1 + lambda*x))/lambda)
}

get.lambda <- function(d){
    median(unlist(apply(d,2,function(tmp){
    bc <- MASS::boxcox(do.call("lm",list(tmp~1)), lambda = seq(.1, .25, .01), plotit = FALSE)
    lam <- bc$x[which.max(bc$y)]
    lam
    })))
}

read.bed <- function (f){
    df=read.table(f, stringsAsFactors=FALSE)
    colnames(df)[1:3]=c('seqnames','start','end')
    makeGRangesFromDataFrame(df[,1:3])
}

read.allbed <- function (conf){
    rbind(read.table(paste0('promoter/',substr(conf$ipReads,1,4)[1],'/row.bed'), stringsAsFactors=FALSE), 
read.table(paste0('enhancer/',substr(conf$ipReads,1,4)[1],'/row.bed'), stringsAsFactors=FALSE))
}

read.counts.trunc <- function (f, trunc=FALSE){
    x = read.counts(f)
    if (trunc)
        x[x==0] = NA
    x
}

read.counts <- function (f){
    if (!file.exists(f)) print(f)
    fh=if(substr(f,nchar(f)-2,nchar(f))=='.gz') gzfile(f) else file(f)
    x=readLines(fh)
    close(fh)
    as.numeric(x)
}

read.data <- function(conf, type){
    list(bed=read.table(paste0(type,'/',substr(conf$ipReads,1,4)[1],'/row.bed'), stringsAsFactors=FALSE), data=do.call("cbind",lapply(paste0(type,'/',conf$ipReads),function(d)read.counts(d))))
}

read.alldata <- function(conf, trunc=FALSE){
    list(bed=rbind(read.table(paste0('promoter/',substr(conf$ipReads,1,4)[1],'/row.bed'), stringsAsFactors=FALSE), 
    read.table(paste0('enhancer/',substr(conf$ipReads,1,4)[1],'/row.bed'), stringsAsFactors=FALSE)),
    data=rbind(do.call("cbind",lapply(paste0('promoter/',conf$ipReads),function(d)read.counts.trunc(d, trunc))), 
    do.call("cbind",lapply(paste0('enhancer/',conf$ipReads),function(d)read.counts.trunc(d, trunc)))))
}

get.rank <- function(k, r){
    x=match(k,r)
    x[!is.na(x)]/length(r)
}

get.generank <- function(r){
    pr=names(sort(r,decreasing=TRUE))
    pr=pr[grepl('_',pr)]
    gsub("_\\d+", "", pr)
}

get.generankscore <- function(r){
    r=sort(r,decreasing=TRUE)
    pr=names(r)
    i=get.uniqgeneid(pr)
    d=r[i]
    names(d)=gsub("_\\d+", "", pr[i])
    d
}

get.rankid <- function(k, r){
    x=match(k,r)
    x=r[x[!is.na(x)]]
    data.frame(gene=x,rank=match(x,r))
}

get.srank <- function(k, r){
    sort(get.rank(k, r))
}

get.gene <- function(gene){
    gene[grepl('_',gene)]
}

get.genename <- function(gene){
    gsub("_\\d+", "", get.gene(gene))
}

get.uniqgeneid <- function(gene){
    match(unique(get.genename(gene)), gsub("_\\d+", "", gene))
}

get.geneid <- function(gene){
    match(get.gene(gene), gene)
}

make.graph <- function(hic, gene, score, weight=1e-20) {
    df=data.frame('Null', setdiff(get.gene(gene),get.gene(hic$gene)),stringsAsFactors=FALSE)
    names(df)=names(hic)
    df=rbind(hic,df)
    g=graph_from_data_frame(df, directed=TRUE)
    j=match(V(g)$name, gene)
    s=score[j]
    s[is.na(s)]=weight
    list(g=g,v=s)
}

get.subgraph <- function(g, gene, uns='_'){
    m=get.edgelist(g)
    induced.subgraph(g,unique(as.character(m[unlist(lapply(paste0(gene,uns), function(d) grep(d,m[,2]))),])))
}

# ported from sSeq
get.auc <- function(yy){
    yy = yy[order(yy)]
    p1 = ecdf(yy)
    1 - trapz(p1(yy), yy)

}

# ported from caTools
trapz <- function (x, y) {
    idx = 2:length(x)
    as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2
}

