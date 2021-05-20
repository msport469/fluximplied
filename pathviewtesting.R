#pathviewtesting
library(pathview)
#BiocManager::install('pathview')
dat<-read.csv('exampledeseqresultdataframe.csv',row.names = c(1))
dat3<-read.csv('exampledeseqresultdataframe.csv')
dat2<-dat[3]


pathview(gene.data = dat[3], 
         pathway.id = '00600',
         species = 'mmu',
         gene.idtype = geneformat,
         kegg.native = F)


lapply(significancetable$keggpathwayid, function(x) pathview(gene.data = inputdat$log2FoldChange, pathway.id = x,
                                                             species = 'mmu',
                                                             gene.idtype = geneformat,
                                                             kegg.native = F)
)
sigkeggpaths<-significancetable$keggpathwayid

pathview(gene.data = inputdat['log2FoldChange'], 
         pathway.id = '00020',
         species = 'mmu',
         gene.idtype = geneformat,
         kegg.native = F)
species
geneformat


lapply(significancetable$keggpathwayid, function(x) pathview(gene.data = inputdat['log2FoldChange'],
                                                             pathway.id = x,
                                                             species = species,
                                                             gene.idtype = geneformat,
                                                             kegg.native = T))
