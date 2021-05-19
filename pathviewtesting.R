#pathviewtesting
library(pathview)
#BiocManager::install('pathview')
dat<-read.csv('exampledeseqresultdataframe.csv',row.names = c(1))
dat3<-read.csv('exampledeseqresultdataframe.csv')
dat2<-dat[3]


pathview(gene.data = dat[3], pathway.id = "00010", species = "mmu", out.suffix = "dat",gene.idtype = 'SYMBOL',kegg.native = F)

str(pv.out)         
pv.out$plot.data.gene$mol.col

data(gse16873.d)


data(rn.list)
names(rn.list)
