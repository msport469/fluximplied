library(org.Mm.eg.db)
library(enrichR)
#dbs <- listEnrichrDbs()
genes<-c('Cpt1a','Idh2','Acadl')
#db<-c('KEGG_2019_Mouse')
res<-enrichr(genes,db)
#kegg<-res[[1]]
RLSgenes<-c('Pfk1','Fbp1','Idh2','Gys1','Pygl','Ppat','Cpsii','Cpsi','ACAC','Cpt1a')
mapIds(org.Mm.eg.db,keys = RLSgenes,column = 'ENTREZID',keytype = 'SYMBOL')
intersect(genes,RLSgenes)





RLSgenes<-c('Pfkm','Pfkl','Fbp1','Idh2','Gys1','Pygl','Ppat','Cad','Cps1','Acaca','Cpt1a')
RLSpathways<-c('glycolysis','glycolysis','gluconeogenesis','TCA cycle','glycogen synthesis','glycogenolysis','de novo purine synthesis',
               'de novo pyrimidine synthesis','urea cycle','fatty acid synthesis','beta oxidation of fatty acids')
RLS<-data.frame(RLSgenes,RLSpathways)
subset<-subset(RLS,RLS$RLSgenes %in% genes)
subset






RLSgenes<-c('Pfkm','Pfkl','Fbp1','Idh2','Gys1','Pygl','Ppat','Cad','Cps1','Acaca','Cpt1a')
RLSpathways<-c('glycolysis','glycolysis','gluconeogenesis','TCA cycle','glycogen synthesis','glycogenolysis','de novo purine synthesis',
               'de novo pyrimidine synthesis','urea cycle','fatty acid synthesis','beta oxidation of fatty acids')
RLS<-data.frame(RLSgenes,RLSpathways)
subset<-subset(RLS,RLS$RLSgenes %in% genes)
colnames(subset)<-c('Rate limiting step genes in your set','Metabolic pathway associated with it')



intersect<-intersect(genes,RLSgenes)
lengthintersect<-length(intersect)
ifzero<-'There are no genes in your set that are in our rate limiting step database. Make sure you used mouse gene symbols like Cpt1a, not the human symbol CPT1A or the ENTREZ ID 1374'
ifnotzero<-print(paste0('Your gene set has --------> ',lengthintersect,' <-------- genes that have been identified as encoding enzymes involved as rate-limiting steps in the gene set you provided'))
ifelse(lengthintersect==0,ifzero,ifnotzero)
subset





RLSgenes<-c('Pfkm','Pfkl','Fbp1','Idh2','Gys1','Pygl','Ppat','Cad','Cps1','Acaca','Cpt1a')
RLSpathways<-c('glycolysis','glycolysis','gluconeogenesis','TCA cycle','glycogen synthesis','glycogenolysis','de novo purine synthesis',
               'de novo pyrimidine synthesis','urea cycle','fatty acid synthesis','beta oxidation of fatty acids')
RLS<-data.frame(RLSgenes,RLSpathways)
subset<-subset(RLS,RLS$RLSgenes %in% genes)
colnames(subset)<-c('Rate limiting step genes in your set','Metabolic pathway associated with it')
intersect<-intersect(genes,RLSgenes)
lengthintersect<-length(intersect)
ifzero<-'There are no genes in your set that are in our rate limiting step database. Make sure you used mouse gene symbols like Cpt1a, not the human symbol CPT1A or the ENTREZ ID 1374'
ifnotzero<-print(paste0('Your gene set has --------> ',lengthintersect,' <-------- genes that have been identified as encoding enzymes involved as rate-limiting steps in the gene set you provided'))
ifelse(lengthintersect==0,ifzero,ifnotzero)
ifelse(lengthintersect==0,'Sorry about that.',subset)



lengthintersect




RLSdatabase<-read.csv('https://raw.githubusercontent.com/msport469/fluximplied/master/RLSdatabase.csv?token=ANC4YV6CMWQ3MBVHDWILXG3ARNO56')




impliedflux <- function(genes) {
  # function to see if there are any rate limiting steps in gene list
  RLSgenes<-c('Pfkm','Pfkl','Fbp1','Idh2','Gys1','Pygl','Ppat','Cad','Cps1','Acaca','Cpt1a')
  RLSpathways<-c('glycolysis','glycolysis','gluconeogenesis','TCA cycle','glycogen synthesis','glycogenolysis','de novo purine synthesis',
                 'de novo pyrimidine synthesis','urea cycle','fatty acid synthesis','beta oxidation of fatty acids')
  RLS<-data.frame(RLSgenes,RLSpathways)
  subset<-subset(RLS,RLS$RLSgenes %in% genes)
  colnames(subset)<-c('Rate limiting step genes in your set','Metabolic pathway associated with it')
  intersect<-intersect(genes,RLSgenes)
  lengthintersect<-length(intersect)
  ifzero<-'There are no genes in your set that are in our rate limiting step database. Make sure you used mouse gene symbols like Cpt1a, not the human symbol CPT1A or the ENTREZ ID 1374'
  ifnotzero<-paste0('Your gene set has --------> ',lengthintersect,' <-------- genes that have been identified as encoding enzymes involved as rate-limiting steps in the gene set you provided')
  
  print(ifelse(lengthintersect==0,ifzero,ifnotzero))
  print(ifelse(lengthintersect==0,'Sorry about that.',subset))
}


impliedflux(c('Idh2','Cpt1a'))

