library(irene)
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
data(markers)
mk=c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3")
nm=c("CLL","Glioma","CRC","B-ALL","CLL-2","MM","PTC")
lx=c(.25,.17,.25,.23,.25,.23,.25)
#conf
conf=read.table('ChIPdesign.txt',sep='\t', header=TRUE)
confs=split(conf, conf$experiment)
#read data
rdata=lapply(1:7,function(i){
read.alldata(confs[[i]])
})
#normalization
ndata=lapply(1:7,function(i){
get.norm.data(confs[[i]],rdata[[i]]$bed,rdata[[i]]$data)
})
#dPCA
res=lapply(1:7,function(i){
meta=confs[[i]]
meta$condition=as.integer(meta$condition=="Healthy")+1
dPCA(meta,rdata[[i]]$bed,rdata[[i]]$data,lambda=lx[i])
})
#ranks
rank=lapply(1:7,function(i){
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$PC[,1]))
page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector
})
#histone marks
marks=lapply(1:7,function(i){
hic=read.table(paste0('in.',substr(confs[[i]]$ipReads,1,4)[1],'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
lapply(1:6,function(j){
g=make.graph(hic[,c('enh','gene')], ndata[[i]]$bed[,4], abs(ndata[[i]]$Dobs[,j]))
page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector
})
})
#nets
nets=lapply(1:7,function(i){
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$PC[,1]))
})
#rewired
rewired=lapply(1:7,function(i){
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
prom=lapply(1:7,function(i){
j=get.geneid(res[[i]]$bed[,4])
gsub("_\\d+", "", res[[i]]$bed[j,4][order(abs(res[[i]]$PC[j,1]),decreasing=TRUE)])
})
save(confs,enh,marks,nets,prom,rank,res,rewired, file="res.rda")
#plot
getAUCs <- function(d) data.frame(x=rep(1:7,each=2),
do.call(rbind,lapply(1:7,function(i) matrix(c(get.auc(get.rank(markers[[i]], get.generank(d[[i]]))), 
get.auc(get.rank(markers$HKG, get.generank(d[[i]]))),
get.auc(get.rank(markers[[i]], prom[[i]])),
get.auc(get.rank(markers$HKG, prom[[i]]))),nrow=2)
)),Type=rep(nm,each=2),Genes=rep(c('Cancer marker genes','Housekeeping genes'),7))
getAUCmat <- function(d) 
do.call(rbind,lapply(1:7,function(i)c(get.auc(get.rank(markers[[i]], get.generank(rank[[i]]))), unlist(lapply(1:6,function(j) get.auc(get.rank(markers[[i]], get.generank(d[[i]][[j]]))))))))
auc=getAUCs(rank)
colnames(auc)[2:3]=c('Irene','Promoter')
dname=function(i) nm[i]
#svglite("f1.svg",6,6)
ggplot(auc, aes(x, color=Genes))+geom_line(aes(y=Irene))+ geom_line(aes(y=Promoter), linetype="dashed")+ scale_x_continuous(breaks=1:7,label=dname)+ labs(x="",y="AUC")+ theme(legend.position=c(.87,.94))
mauc=getAUCmat(marks)
colnames(mauc)=c('Irene',mk)


