# Irene: Integrative Ranking with Epigenetic Network of Enhancers

## Installation
```r
library(devtools)
devtools::install_github("uni-hd/Irene")
```

## Quick start 
Rerun all the test cases in the paper (in R session): 
```r
library(irene)
download.file("https://github.com/uni-hd/Irene-data/archive/master.zip","Irene-data-master.zip")
unzip("Irene-data-master.zip")
setwd("Irene-data-master")
source(system.file("extdata", "example.R", package="irene"), echo=TRUE)
```
Output the Irene rank list for the first test case
```r
write.table(format(get.generankscore(rank[[1]]), digits=3),file="PR.txt",sep="\t",quote=F,col.names=F)
```

## Step-by-step instructions
### Input files
A four-columns meta file (ChIPdesign.txt in our example) with the following header: experiment, factor, condition, ipReads.
* experiment: the study ID, should be the same for the healthy and diseased tissues. 
* factor: the histone marks. 
* condition: healthy and diseased states.
* ipReads: the file which contains the sum of ChIP-Seq signals over the pre-defined promoter ([hg19](https://github.com/uni-hd/Irene-data/blob/master/promoter/hg19/row.bed), [hg38](https://github.com/uni-hd/Irene-data/blob/master/promoter/hg38/row.bed)) and enhancer ([hg19](https://github.com/uni-hd/Irene-data/blob/master/enhancer/hg19/row.bed),[hg38](https://github.com/uni-hd/Irene-data/blob/master/enhancer/hg38/row.bed)) regions. The numbers can be obtained using [bigWigAverageOverBed](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/). 

A promoter-enhancer interaction table with the following header: chromosome names, start, end, gene names, enhancer names. It can be computed from Hi-C and/or ChIA-PET experiments. We have prepared the [hg19](https://github.com/uni-hd/Irene-data/blob/master/in.hg19.bed.gz) and [hg38](https://github.com/uni-hd/Irene-data/blob/master/in.hg38.bed.gz) files from the [4DGenome database](https://4dgenome.research.chop.edu/Download.html).
### R settings
```r
library(irene)
Sys.setlocale("LC_NUMERIC","C")
options(stringsAsFactors = FALSE)
```
### Load input data
```r
conf=read.table('ChIPdesign.txt',sep='\t', header=TRUE)
meta=conf[conf$experiment==1,]
meta$condition=as.integer(meta$condition=="Healthy")+1
rdata=read.alldata(meta)
```
### Compute dPCA
```r
res=dPCA(meta,rdata$bed,rdata$data,lambda=0.25)
```
### Compute gene ranks
```r
hic=read.table('in.hg38.bed.gz',sep='\t',col.names=c('seqnames','start','end','gene','enh'))
g=make.graph(hic[,c('enh','gene')], res$bed[,4], abs(res$PC[,1]))
rank=page_rank(g$g, algo="arpack", personalized=g$v, directed=TRUE)$vector
write.table(format(get.generankscore(rank), digits=3),file="PR.txt",sep="\t",quote=F,col.names=F)
```
### Promoter-enhancer network permutations
You can also replace the promoter-enhancer interactions with nearest enhancer assignment, 
```r
hic=read.table('nearest.hg38.bed.gz',sep='\t',col.names=c('seqnames','start','end','gene','enh'))
```
or create a randomly shuffled network to see what the ranks may look like. 
```r
g$g=rewire(g$g, with=each_edge(prob=1))
```
### Get the rank lists
```r
j=get.geneid(res$bed[,4])
prom=gsub("_\\d+", "", res$bed[j,4][order(abs(res$PC[j,1]),decreasing=TRUE)])
enh=get.generank(rank)
```
### Export reports as a local website
Put "index.Rmd" and "\_site.yml" in the current working directory, then run render_site from the _rmarkdown_ package
```r
res=list(res)
hic=list(hic)
enh=list(enh)
prom=list(prom)
rank=list(rank)
rmarkdown::render_site()
```
It creates an HTML file which you can open with Internet browsers. 

### Also check out our results:
[https://uni-hd.github.io/Irene/](https://uni-hd.github.io/Irene/)
