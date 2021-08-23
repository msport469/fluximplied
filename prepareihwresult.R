#load in deseqresult from previous publication (Frontiers in Immunology 2021), then we'll save it
#saveRDS(ihw.res.DPvsDN.all,'DPDNdfmouse.RDS')
inputdat<-readRDS('DPDNdfmouse.RDS')
inputdatgsea<-inputdat
inputdat$gene<-rownames(inputdat)
#subset on p values
inputdatgsea <- inputdatgsea[inputdatgsea[,"weighted_pvalue"]<0.05,]
#differentiate up and down
inputdatgseaup <- inputdatgsea[inputdatgsea[,"log2FoldChange"]>0.5,]
inputdatgseadown <- inputdatgsea[inputdatgsea[,"log2FoldChange"]<(-0.5),]
#save genes as vectors.
upgenes <- rownames(inputdatgseaup)
downgenes <- rownames(inputdatgseadown)