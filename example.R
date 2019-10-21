library(irene)
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)


#conf
confs=lapply(1:9,function(i){
conf=read.table('ChIPdesign.txt',sep='\t',stringsAsFactors=FALSE,col.names=c('experiment','factor','condition','ipReads','ctReads','peaks'))
conf=conf[conf$experiment==i,]
conf[conf$factor!="H3K9ac",]
})
#dPCA
res=lapply(1:9,function(i){
data=read.alldata(confs[[i]])
dPCA(confs[[i]],data$bed,data$data,verbose=FALSE)
#par(mfrow=c(2,3))
#for(i in 1:6) barplot(prcomp(res$Dobs)$rotation[,i],names=mk,main = paste0("PC",i),las=2)
})
#H3K36me3
order2 <- function(df) df[order(abs(df[,ncol(df)]),decreasing=TRUE),]
k36=lapply(1:9,function(i){
data=read.data(confs[[i]],'promoter')
order2(cbind(data$bed, dPCA(confs[[i]],data$bed,data$data,verbose=FALSE)$PC))
})
#bed
bed=lapply(1:9,function(i){
read.allbed(confs[[i]])
})
#ranks
rank=lapply(1:9,function(i){
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$PC[,1]))
#g$g=rewire(g$g, with=each_edge(prob=1))
page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)
})
#rewired
rewired=lapply(1:9,function(i){
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$PC[,1]))
lapply(1:100,function(j){
g$g=rewire(g$g, with=each_edge(prob=1))
get.generank(page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE))
})
})
#gene ranks
enh=sapply(rank,get.generank)
prom=lapply(1:9,function(i){
j=get.geneid(res[[i]]$bed[,4])
gsub("_\\d+", "", res[[i]]$bed[j,4][order(abs(res[[i]]$PC[j,1]),decreasing=TRUE)])
})
#ECDF AUC
auc=do.call("rbind",lapply(1:9,function(i){
require(sSeq)
png(tempfile(paste0(i,".png")))
x=ecdfAUC(data.frame(Irene=get.srank(markers[[i]], enh[[i]]), Promoter=get.srank(markers[[i]], prom[[i]]), K36=get.srank(markers[[i]], get.genename(k36[[i]][,4]))),lineType=c(1,1,1),lwd=2,xlab="Rank (%)")
dev.off()
x
}))
summary(auc[,1]-auc[,2])

