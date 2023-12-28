#this script exists to update the RLS database
#load libraries
library(org.Mm.eg.db)
library(org.Hs.eg.db)
#load our already published database
RLSdatabase<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplieddev/master/RLSdatabase.csv',
                      stringsAsFactors = F)
#look at the top
head(RLSdatabase)
#if necessary, write the csv and manually edit/curate it to add new rate limiting steps.
#map the IDs to create a column of mouse ENTREZIDs and human ENTREZIDs
RLSdatabase$mouse.entrez <- mapIds(org.Mm.eg.db,keys = RLSdatabase$mouse.gene.symbol,
                                   keytype = 'SYMBOL',
                                   column = 'ENTREZID')
RLSdatabase$human.entrez <- mapIds(org.Hs.eg.db,keys = RLSdatabase$human.gene.symbol,
                                   keytype = 'SYMBOL',
                                   column = 'ENTREZID')
#look at the top
head(RLSdatabase)
#rewrite the csv so I can push it to github later and publish it.
write.csv(RLSdatabase,'RLSdatabase.csv')

sessionInfo()