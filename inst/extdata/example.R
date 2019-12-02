library(irene)
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
data(markers)
mk=c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3")
nm=c("CLL","Glioma","CRC","B-ALL","mCLL","MM","PTC")
lx=c(.25,.17,.25,.23,.25,.23,.25)
#conf
conf=read.table('ChIPdesign.txt',sep='\t', header=TRUE)
confs=split(conf, conf$experiment)
#read data
print('reading data ...')
rdata=lapply(1:7,function(i){
read.alldata(confs[[i]], trunc=T)
})
#normalization
print('normalization ...')
ndata=lapply(1:7,function(i){
get.norm.data(confs[[i]],rdata[[i]]$bed,rdata[[i]]$data,lambda=lx[i],transform=T)
})
#dPCA
print('computing dPCA ...')
res=lapply(1:7,function(i){
meta=confs[[i]]
meta$condition=as.integer(meta$condition=="Healthy")+1
dPCA(meta,rdata[[i]]$bed,rdata[[i]]$data,lambda=lx[i])
})
#ranks
print('computing ranks 1/4 ...')
rank=lapply(1:7,function(i){
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$PC[,1]))
page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector
})
#histone marks
print('computing ranks 2/4 ...')
marks=lapply(1:7,function(i){
hic=read.table(paste0('in.',substr(confs[[i]]$ipReads,1,4)[1],'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
lapply(1:6,function(j){
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$Dobs[,j]))
page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector
})
})
#nearest enhancer
print('computing ranks 3/4 ...')
nearest=lapply(1:7,function(i){
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('nearest.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$PC[,1]))
page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector
})
#rewired
print('computing ranks 4/4 ...')
rewired=lapply(1:7,function(i){
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(res[[i]]$PC[,1]))
lapply(1:10,function(j){
g$g=rewire(g$g, with=each_edge(prob=1))
get.generank(page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector)
})
})
#gene ranks
enh=lapply(rank,get.generank)
promrank=lapply(1:7,function(i){
j=get.geneid(res[[i]]$bed[,4])
nv=abs(res[[i]]$PC[j,1])
names(nv)=res[[i]]$bed[j,4]
nv
})
prom=lapply(1:7,function(i){
j=get.geneid(res[[i]]$bed[,4])
gsub("_\\d+", "", res[[i]]$bed[j,4][order(abs(res[[i]]$PC[j,1]),decreasing=TRUE)])
})
#plot
print('plotting figures ...')
library(ggplot2)
library(reshape2)
dframe=function(d,col.names=NULL,row.names=NULL,m.colnames=NULL,m.rownames=NULL,melt=FALSE){
if(!is.null(col.names)) colnames(d)=col.names
if(!is.null(row.names)) rownames(d)=row.names
if(melt) d=melt(d)
if(!is.null(m.colnames)) colnames(d)=m.colnames
if(!is.null(m.rownames)) rownames(d)=m.rownames
d
}
v=c(2,5,6,4,1,7,3)
hic=lapply(1:7,function(i) read.table(paste0('in.',substr(confs[[i]]$ipReads,1,4)[1],'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','prom','enh')))
dname=function(i) nm[v][i]
ad <- function(x,y,z,d) data.frame(Type=x,Method=y,Genes=z,value=sort(d))
auc=dframe(do.call(rbind,lapply(v,function(i)c(get.auc(get.rank(markers[[i]], get.generank(rank[[i]]))), unlist(lapply(1:6,function(j) get.auc(get.rank(markers[[i]], get.generank(marks[[i]][[j]])))))))),col.names=c('dPC1',mk),m.colnames=c('Type','Mark','value'),melt=T)
ggplot(auc, aes(x=Type,y=value,colour=Mark,shape=Mark))+geom_point(size=3)+scale_shape_manual(values=c(19,17,12:8))+scale_color_manual(values=RColorBrewer::brewer.pal(7,"Dark2"))+ scale_x_continuous(breaks=1:7,label=dname)+ labs(x="",y="AUC")
ggsave("f1.pdf",width=6,height=6,units="in")
df=do.call(rbind,lapply(1:7,function(i) rbind(
ad(nm[i],'Irene','Cancer marker genes',get.rank(markers[[i]], get.generank(rank[[i]]))), 
ad(nm[i],'Irene','Housekeeping genes',get.rank(markers$HKG, get.generank(rank[[i]]))),
ad(nm[i],'Promoter','Cancer marker genes',get.rank(markers[[i]], prom[[i]])),
ad(nm[i],'Promoter','Housekeeping genes',get.rank(markers$HKG, prom[[i]])))))
ggplot(df,aes(value, colour=Method, linetype=Genes))+stat_ecdf(pad=TRUE)+labs(x="Rank", y="ECDF") +facet_wrap(.~Type)+ theme(legend.position=c(0.8,0.14))
ggsave("f2.pdf",width=6,height=6,units="in")
df=do.call(rbind,lapply(1:7,function(i) rbind(
ad(nm[i],'Promoter','Cancer marker genes',get.generankscore(promrank[[i]])[markers[[i]]]),
ad(nm[i],'Promoter','Housekeeping genes',get.generankscore(promrank[[i]])[markers$HKG]))))
ggplot()+geom_boxplot(data=df, aes(Type, value, colour=Genes),outlier.shape = NA) +labs(x='', y="dPC1")+ theme(legend.position=c(.85,.9))
ggsave("f6.pdf",width=6,height=6,units="in")
df=dframe(do.call(rbind,lapply(1:7,function(i)summary(prcomp(res[[i]]$Dobs))$importance[2,])),row.names=nm,col.names=paste0('dPC',1:6),m.colnames=c('Type','dPC','value'),melt=T)
ggplot(df,aes(Type, value, fill=dPC))+geom_bar(stat="identity",position="stack")+labs(x='', y="Variance")+scale_y_continuous(labels = scales::percent) 
ggsave("f3.pdf",width=6,height=6,units="in")
df = do.call(rbind,lapply(1:7, function(i) data.frame(melt(prcomp(res[[i]]$Dobs)$rotation),rep(mk,6),nm[i])))
colnames(df)=c("id","PC","Loadings","Mark","type")
ggplot(df, aes(x=Mark, y=Loadings, fill=Mark)) + geom_bar(stat='identity')+facet_grid(PC~type) +coord_flip()+theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("f4.pdf",width=6,height=6,units="in")
df=dframe(do.call(rbind,lapply(v,function(i) unlist(lapply(1:10,function(j) get.auc(get.rank(markers[[i]], rewired[[i]][[j]])))))),row.names=nm[v],melt=T)
au=dframe(do.call(rbind,lapply(v,function(i)c(get.auc(get.rank(markers[[i]], get.generank(rank[[i]]))),get.auc(get.rank(markers[[i]], prom[[i]])),get.auc(get.rank(markers[[i]], get.generank(nearest[[i]])))))),col.names=c('Irene','Promoter','Nearest'),m.colnames=c('Type','Method','value'),melt=T)
ggplot() + geom_boxplot(data=df, aes(Var1, value),outlier.shape = NA) + geom_line(data=au, aes(Type, value, colour=Method))+ labs(x='',y="AUC")+ theme(legend.position=c(.85,.89))
ggsave("f5.pdf",width=6,height=6,units="in")

