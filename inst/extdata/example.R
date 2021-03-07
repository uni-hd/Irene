library(irene)
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
data(markers)
mk=c("H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3")
nm=c("CLL","Glioma","CRC","B-ALL","mCLL","MM","PTC")
lx=rep(.25,7)
set.seed(1)
#conf
conf=read.table('ChIPdesign.txt',sep='\t', header=TRUE)
confs=split(conf, conf$experiment)
#read data
print('reading data ...')
rdata=lapply(1:7,function(i){
d = read.alldata(confs[[i]], trunc=TRUE)
list(bed=d$bed, data=d$data)
})
#dPCA
print('computing dPCA ...')
res=lapply(1:7,function(i){
meta=confs[[i]]
meta$condition=as.integer(meta$condition=="Healthy")+1
dPCA(meta,rdata[[i]]$bed,rdata[[i]]$data,lambda=lx[i])
})
#merge PC1 and PC2
for(i in 1:7) res[[i]]$PC[,1]=abs(res[[i]]$PC[,1])+abs(res[[i]]$PC[,2])
#precompiled list for tissue-specific enhancers
ji= lapply(1:7,function(i){
id= read.table(paste0(i,'.id'))[,1]
d = res[[i]]$bed[,4]
grepl('_\\d+',d) | (d %in% id)
})
#ranks
print('computing ranks 1/3 ...')
rank=lapply(1:7,function(i){j=ji[[i]]
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[j,4], abs(res[[i]]$PC[j,1]))
page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector
})
#nearest enhancer
print('computing ranks 2/3 ...')
nearest=lapply(1:7,function(i){j=ji[[i]]
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('nearest.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[j,4], abs(res[[i]]$PC[j,1]))
page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector
})
#rewired iterations
ri = 1:3
print('computing ranks 3/3 ...')
rewired=lapply(1:7,function(i){j=ji[[i]]
spec=substr(confs[[i]]$ipReads,1,4)[1]
hic=read.table(paste0('in.',spec,'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res[[i]]$bed[j,4], abs(res[[i]]$PC[j,1]))
lapply(ri,function(j){
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
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(2)
shortnames <- function(d){
d$CRC=d$`Colorectal cancer`
d$PTC=d$`Papillary thyroid carcinoma`
d$MM=d$`Multiple myeloma`
d$mCLL=d$CLL
d}
#tissue-specific genes
ts=c('COLON (BULK TISSUE)','THYROID (BULK TISSUE)','CD19+ B CELLS','BRAIN (BULK)')
download.file("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=ARCHS4_Tissues","ARCHS4_Tissues.txt")
x=read.table("ARCHS4_Tissues.txt",sep='\t',row.names=1)
x=split(x[,-ncol(x)],rownames(x))
tissues=x[c(ts[c(3,4,1,3,3,3,2)])]
names(tissues)=names(markers)[1:7]
markers = shortnames(markers)
tissues = shortnames(tissues)
dname=function(i) nm[v][i]
ad <- function(x,y,z,d) data.frame(Type=x,Method=y,Genes=z,value=sort(d))
df=do.call(rbind,lapply(1:7,function(i) rbind(
ad(nm[i],'Irene',   'Cancer marker genes',  get.rank(markers[[i]], enh[[i]])), 
ad(nm[i],'Irene',   'Tissue-specific genes',get.rank(tissues[[i]], enh[[i]])),
ad(nm[i],'Promoter','Cancer marker genes',  get.rank(markers[[i]], prom[[i]])),
ad(nm[i],'Promoter','Tissue-specific genes',get.rank(tissues[[i]], prom[[i]])))))
au1=data.frame(label=round(unlist(lapply(1:7,function(i) get.auc(get.rank(markers[[nm[i]]], enh[[i]])))),digits=2),Type=nm, Method='Irene',Genes='Cancer marker genes')
au2=data.frame(label=round(unlist(lapply(1:7,function(i) get.auc(get.rank(markers[[nm[i]]], prom[[i]])))),digits=2),Type=nm, Method='Promoter',Genes='Cancer marker genes')
au3=data.frame(label=round(unlist(lapply(1:7,function(i) get.auc(get.rank(tissues[[nm[i]]], enh[[i]])))),digits=2),Type=nm, Method='Irene',Genes='Tissue-specific genes')
au4=data.frame(label=round(unlist(lapply(1:7,function(i) get.auc(get.rank(tissues[[nm[i]]], prom[[i]])))),digits=2),Type=nm, Method='Promoter',Genes='Tissue-specific genes')
p=ggplot(df,aes(value, colour=Method, linetype=Genes))+stat_ecdf(pad=TRUE)+labs(x="Rank", y="ECDF") +facet_wrap(.~Type)+ theme(legend.position=c(0.8,0.14))
p=p+geom_text(data.frame(label='AUCs',Type=nm, Method='Irene',Genes='Cancer marker genes'), mapping = aes(x = 0.92, y = 0.33, label = label), size=3, show.legend=FALSE)
p=p+geom_text(au1, mapping = aes(x = 0.92, y = 0.25, label = label), size=2.5, show.legend=FALSE)+geom_segment(aes(x = 0.7, y = 0.25, xend = 0.85, yend = 0.25), linetype="solid", colour = cols[1])
p=p+geom_text(au2, mapping = aes(x = 0.92, y = 0.20, label = label), size=2.5, show.legend=FALSE)+geom_segment(aes(x = 0.7, y = 0.20, xend = 0.85, yend = 0.20), linetype="solid", colour = cols[2])
p=p+geom_text(au3, mapping = aes(x = 0.92, y = 0.15, label = label), size=2.5, show.legend=FALSE)+geom_segment(aes(x = 0.7, y = 0.15, xend = 0.85, yend = 0.15), linetype="dashed", colour = cols[1])
p=p+geom_text(au4, mapping = aes(x = 0.92, y = 0.10, label = label), size=2.5, show.legend=FALSE)+geom_segment(aes(x = 0.7, y = 0.10, xend = 0.85, yend = 0.10), linetype="dashed", colour = cols[2])
ggsave(p,file="f2.pdf",width=6,height=6,units="in")
df=do.call(rbind,lapply(1:7,function(i) {
hic=read.table(paste0('in.',substr(confs[[i]]$ipReads,1,4)[1],'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','prom','enh'))
rbind(
ad(nm[i],'Promoter','Cancer marker genes',  abs(res[[i]]$PC[gsub("_\\d+", "", res[[i]]$bed[,4]) %in% markers[[i]],1])),
ad(nm[i],'Promoter','Tissue-specific genes',abs(res[[i]]$PC[gsub("_\\d+", "", res[[i]]$bed[,4]) %in% tissues[[i]],1])))}))
p=ggplot()+geom_boxplot(data=df, aes(Type, value, colour=Genes),outlier.shape = NA) +labs(x='', y="dPC1")+ theme(legend.position=c(.85,.9))
df=rbind(df,do.call(rbind,lapply(1:7,function(i) {
hic=read.table(paste0('in.',substr(confs[[i]]$ipReads,1,4)[1],'.bed.gz'),sep='\t',col.names=c('seqnames','start','end','prom','enh'))
rbind(
ad(nm[i],'Enhancer','Cancer marker genes',  abs(res[[i]]$PC[res[[i]]$bed[,4] %in% hic[get.genename(hic$prom) %in% markers[[i]],5],1])),
ad(nm[i],'Enhancer','Tissue-specific genes',abs(res[[i]]$PC[res[[i]]$bed[,4] %in% hic[get.genename(hic$prom) %in% tissues[[i]],5],1])))})))
p=ggplot()+facet_grid(Method~.,scales="free")+geom_boxplot(data=df, aes(Type, value, colour=Genes),outlier.shape = NA) +labs(x='', y="dPC1")+ theme(legend.position=c(.8,.9))
ggsave("f1.pdf",width=6,height=6,units="in")
df=dframe(do.call(rbind,lapply(1:7,function(i)summary(prcomp(res[[i]]$Dobs))$importance[2,])),row.names=nm,col.names=paste0('dPC',1:6),m.colnames=c('Type','dPC','value'),melt=T)
p=ggplot(df,aes(Type, value, fill=dPC))+geom_bar(stat="identity",position="stack")+labs(x='', y="Variance")+scale_y_continuous(labels = scales::percent) 
ggsave("f3.pdf",width=6,height=6,units="in")
df = do.call(rbind,lapply(1:7, function(i) data.frame(melt(prcomp(res[[i]]$Dobs)$rotation),rep(mk,6),nm[i])))
colnames(df)=c("id","PC","Loadings","Mark","type")
p=ggplot(df, aes(x=Mark, y=Loadings, fill=Mark)) + geom_bar(stat='identity')+facet_grid(PC~type) +coord_flip()+theme(legend.position="none",axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("f4.pdf",width=6,height=6,units="in")
df=dframe(do.call(rbind,lapply(v,function(i) unlist(lapply(ri,function(j) get.auc(get.rank(markers[[i]], rewired[[i]][[j]])))))),row.names=nm[v],melt=T)
au=dframe(do.call(rbind,lapply(v,function(i)c(get.auc(get.rank(markers[[i]], get.generank(rank[[i]]))),get.auc(get.rank(markers[[i]], prom[[i]])),get.auc(get.rank(markers[[i]], get.generank(nearest[[i]])))))),col.names=c('Irene','Promoter','Nearest'),m.colnames=c('Type','Method','value'),melt=T)
p=ggplot(df, aes(Var1, value)) + geom_violin() +geom_boxplot(width=0.2,coef=0,outlier.shape=NA) + geom_line(data=au, aes(Type, value, colour=Method))+ labs(x='',y="AUC")+ theme(legend.position=c(.85,.89))
ggsave("f5.pdf",width=6,height=6,units="in")
#validation of AUCs on multiple cancer genesets
info=data.frame(file=c('cancermine_collated.tsv.gz','Compendium_Cancer_Genes.tsv.gz'),
	gene=c('V7','V1'),type=c('V4','V4'),gset=c('CancerMine','IntOGen'),out=c('cancermine','intogen'),size=c(20,20))
for(k in 1:2){d=info[k,]
x=read.table(gzfile(d$file),sep='\t',skip=1)
x=unique(x)
x=split(x,x[,d$type])
x=x[unlist(lapply(x,nrow))>d$size]
df=do.call(rbind,lapply(1:7,function(i){
gr=get.generank(rank[[i]])
set.seed(1)
d1=do.call(rbind,lapply(1:length(x),function(j) c(nrow(x[[j]]),get.auc(get.rank(x[[j]][,d$gene], gr)))))
d2=do.call(rbind,lapply(1:length(x),function(j) c(nrow(x[[j]]),get.auc(get.rank(sample(setdiff(tissues[[i]],x[[j]][,d$gene]),nrow(x[[j]]),replace=F), gr)))))
df=data.frame(d1,GeneSet=d$gset)
df=rbind(df,data.frame(d2,GeneSet='Tissue-specific genes'))
data.frame(df,Type=nm[i])
}))
df$GeneSet=factor(df$GeneSet,levels=c(d$gset,'Tissue-specific genes'))
p=ggplot(df,aes(x=GeneSet, y=X2, colour=GeneSet)) + geom_violin() +geom_boxplot(width=0.2,coef=0,outlier.shape=NA)+  facet_grid(.~Type)+labs(x='Test cases', y="AUCs")+ theme(axis.text.x = element_blank(), legend.position=c(.86,.07), legend.text=element_text(size=8),legend.title = element_blank(),legend.margin=margin(t=0, unit='cm'))+ scale_size(guide="none")
ggsave(p,file=paste0(d$out,'.pdf'),height=6,width=6)
}
