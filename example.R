library(irene)
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
mk=c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3")

#conf
confs=lapply(1:9,function(i){
conf=read.table('ChIPdesign.txt',sep='\t',stringsAsFactors=FALSE,col.names=c('experiment','factor','condition','ipReads','ctReads','peaks'))
conf=conf[conf$experiment==i,]
conf[conf$factor!="H3K9ac",]
})
#normalization
res=lapply(1:9,function(i){
data=read.alldata(confs[[i]])
get.norm.data(confs[[i]],data$bed,data$data)
})
#dPCA
dres=lapply(1:9,function(i){
data=read.alldata(confs[[i]])
dPCA(confs[[i]],data$bed,data$data)
})
#histone marks
marks=lapply(1:9,function(i){
hic=read.table(paste0('in.',substr(confs[[i]]$ipReads,1,4)[1],'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
lapply(1:6,function(j){
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$Dobs[,j]))
page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector
})
})
#nets
nets=lapply(1:9,function(i){
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$PC[,1]))
})
#ranks
rank=lapply(1:9,function(i){
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$PC[,1]))
page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector
})
#rewired
rewired=lapply(1:9,function(i){
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$PC[,1]))
lapply(1:100,function(j){
g$g=rewire(g$g, with=each_edge(prob=1))
get.generank(page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector)
})
})
#gene ranks
enh=sapply(rank,get.generank)
prom=lapply(1:9,function(i){
j=get.geneid(res[[i]]$bed[,4])
gsub("_\\d+", "", res[[i]]$bed[j,4][order(abs(res[[i]]$Dobs[j,1]),decreasing=TRUE)])
})
save(confs,enh,marks,nets,prom,rank,res,rewired, file="res.rda")


