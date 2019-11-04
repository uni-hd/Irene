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
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[,4], abs(boxCox(res[[i]]$Dobs[,j],lx[i])))
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
#test plot
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
getECDF <- function(d,i) do.call(rbind,lapply(1:6,function(j) as.matrix(data.frame(nm[i],mk[j],sort(get.rank(markers[[i]], get.generank(d[[j]])))))))
getECDFmat <- function() 
df=rbind(do.call(rbind,lapply(1:7,function(i)getECDF(marks[[i]],i))),
do.call(rbind,lapply(1:7,function(i)as.matrix(data.frame(nm[i],'Irene',sort(get.rank(markers[[i]], get.generank(rank[[i]]))))))))
df=as.data.frame(getECDFmat())
colnames(df)=c('Type','Method','value')
df$Type=as.factor(df$Type)
df$Method=as.factor(df$Method)
df$value=as.numeric(df$value)
#svglite("f2.svg",6,6)
ggplot(df,aes(value, colour=Method))+stat_ecdf(pad = FALSE)+labs(x="Rank", y="ECDF") +facet_wrap(.~Type)+ theme(legend.position=c(0.8,0.15))
#plot
dframe=function(d,col.names=NULL,row.names=NULL,m.colnames=NULL,m.rownames=NULL,melt=FALSE){
if(!is.null(col.names)) colnames(d)=col.names
if(!is.null(row.names)) rownames(d)=row.names
if(melt) d=melt(d)
if(!is.null(m.colnames)) colnames(d)=m.colnames
if(!is.null(m.rownames)) rownames(d)=m.rownames
d
}
v=c(2,5,6,4,1,7,3)
dname=function(i) nm[v][i]
ad <- function(x,y,z,d) data.frame(Type=x,Method=y,Genes=z,value=sort(d))
#svglite("f1.svg",6,6)
auc=dframe(do.call(rbind,lapply(v,function(i)c(get.auc(get.rank(markers[[i]], get.generank(rank[[i]]))), unlist(lapply(1:6,function(j) get.auc(get.rank(markers[[i]], get.generank(marks[[i]][[j]])))))))),col.names=c('dPC1',mk),m.colnames=c('Type','Mark','value'),melt=T)
ggplot(auc, aes(x=Type,y=value,colour=Mark,shape=Mark))+geom_point(size=3)+scale_shape_manual(values=c(19,17,12:8))+scale_color_manual(values=RColorBrewer::brewer.pal(7,"Dark2"))+ scale_x_continuous(breaks=1:7,label=dname)+ labs(x="",y="AUC")
#svglite("f2.svg",6,6)
df=do.call(rbind,lapply(1:7,function(i) rbind(
ad(nm[i],'Irene','Cancer marker genes',get.rank(markers[[i]], get.generank(rank[[i]]))), 
ad(nm[i],'Irene','Housekeeping genes',get.rank(markers$HKG, get.generank(rank[[i]]))),
ad(nm[i],'Promoter','Cancer marker genes',get.rank(markers[[i]], prom[[i]])),
ad(nm[i],'Promoter','Housekeeping genes',get.rank(markers$HKG, prom[[i]])))))
ggplot(df,aes(value, colour=Method, linetype=Genes))+stat_ecdf(pad=TRUE)+labs(x="Rank", y="ECDF") +facet_wrap(.~Type)+ theme(legend.position=c(0.8,0.14))
#svglite("f3.svg",6,6)
df=dframe(do.call(rbind,lapply(1:7,function(i)res[[i]]$proj[2,])),row.names=nm,col.names=paste0('dPC',1:6),m.colnames=c('Type','dPC','value'),melt=T)
ggplot(df,aes(Type, value, fill=dPC))+geom_bar(stat="identity",position="stack")+labs(x='', y="Variance")+scale_y_continuous(labels = scales::percent) 
#svglite("f4.svg",6,6)
df = do.call(rbind,lapply(1:7, function(i) data.frame(melt(prcomp(res[[i]]$Dobs)$rotation),rep(mk,6),nm[i])))
colnames(df)=c("id","PC","Loadings","Mark","type")
ggplot(df, aes(x=Mark, y=Loadings, fill=Mark)) + geom_bar(stat='identity')+facet_grid(PC~type) +coord_flip()+theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5))
#svglite("f5.svg",6,6)
df=dframe(do.call(rbind,lapply(v,function(i) unlist(lapply(1:10,function(j) get.auc(get.rank(markers[[i]], rewired[[i]][[j]])))))),row.names=nm[v],melt=T)
au=dframe(do.call(rbind,lapply(v,function(i)c(get.auc(get.rank(markers[[i]], get.generank(rank[[i]]))),get.auc(get.rank(markers[[i]], prom[[i]]))))),col.names=c('Irene','Promoter'),m.colnames=c('Type','Method','value'),melt=T)
ggplot() + geom_boxplot(data=df, aes(Var1, value),outlier.shape = NA) + geom_line(data=au, aes(Type, value, linetype=Method),col="red")+ labs(x='',y="AUC")

bed=res[[i]]$bed
rownames(bed)=res[[i]]$bed[,4]
hi=cbind(hic,bed[hic$gene,])
hi=hi[hi$seqnames==hi$V1,]
hlen=data.frame(hi[,c('gene','enh')],hi[,c('start','end')]-rowMeans(hi[,c('V2','V3')]))
hlen$gene=get.genename(hlen$gene)


