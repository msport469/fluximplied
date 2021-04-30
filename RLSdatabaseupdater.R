#this script exists to update the RLS database
library(org.Mm.eg.db)
library(org.Hs.eg.db)
RLSdatabase<-read.csv('https://raw.githubusercontent.com/msport469/fluximplied/master/RLSdatabase.csv?token=ANC4YV2BBAO6ACQMKN6WLFLARNP7G',
                      stringsAsFactors = F)
head(RLSdatabase)

RLSdatabase$mouse.entrez <- mapIds(org.Mm.eg.db,keys = RLSdatabase$mouse.gene.symbol,
                                   keytype = 'SYMBOL',
                                   column = 'ENTREZID')
RLSdatabase$human.entrez <- mapIds(org.Hs.eg.db,keys = RLSdatabase$human.gene.symbol,
                                   keytype = 'SYMBOL',
                                   column = 'ENTREZID')
write.csv(RLSdatabase,'RLSdatabase.csv')
