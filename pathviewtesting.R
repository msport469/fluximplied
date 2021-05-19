#pathviewtesting
library(pathview)
#BiocManager::install('pathview')
dat<-read.csv('exampledeseqresultdataframe.csv',row.names = c(1))
dat3<-read.csv('exampledeseqresultdataframe.csv')
dat2<-dat[3]


pathview(gene.data = dat[3], pathway.id = "00010",
         species = species,
         gene.idtype = geneformat,
         kegg.native = F)
