library(ggplot2)
library(viridis)
library(enrichR)
library(gridExtra)
#load in deseqresult
saveRDS(ihw.res.DPvsDN.all,'DPDNdf.RDS')
inputdat<-readRDS('DPDNdf.RDS')
inputdatgsea<-inputdat
inputdat$gene<-rownames(inputdat)
inputdatgsea <- inputdatgsea[inputdatgsea[,"weighted_pvalue"]<0.05,]
inputdatgseaup <- inputdatgsea[inputdatgsea[,"log2FoldChange"]>0,]
inputdatgseadown <- inputdatgsea[inputdatgsea[,"log2FoldChange"]<0,]


upgenes <- rownames(inputdatgseaup)
downgenes <- rownames(inputdatgseadown)
