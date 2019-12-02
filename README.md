# Irene: Integrative Ranking with Epigenetic Network of Enhancers

## Usage (in R session): 
```r
library(devtools)
devtools::install_github("uni-hd/Irene")
download.file("https://github.com/uni-hd/Irene-data/archive/master.zip","Irene-data-master.zip")
unzip("Irene-data-master.zip")
setwd("Irene-data-master")
source(system.file("extdata", "example.R", package="irene"), echo=TRUE)

#output Irene rank list for the first test case
write.table(format(get.generankscore(rank[[1]]), digits=3),file="PR.txt",sep="\t",quote=F,col.names=F)
```

## Also check out our results:
[https://uni-hd.github.io/Irene/](https://uni-hd.github.io/Irene/)
