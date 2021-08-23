#load in deseqresult from previous publication (Frontiers in Immunology 2021), then we'll save it
#saveRDS(ihw.res.DPvsDN.all,'DPDNdfmouse.RDS')
inputdatgseamouse<-readRDS('DPDNdfmouse.RDS')
inputdatgseamouse$gene<-rownames(inputdatgseamouse)
#subset on p values
inputdatgseamouse <- inputdatgseamouse[inputdatgseamouse[,"weighted_pvalue"]<0.05,]
#differentiate up and down
inputdatgseamouseup <- inputdatgseamouse[inputdatgseamouse[,"log2FoldChange"]>0.5,]
inputdatgseamousedown <- inputdatgseamouse[inputdatgseamouse[,"log2FoldChange"]<(-0.5),]
#save genes as vectors.
upgenesmouse <- rownames(inputdatgseamouseup)
downgenesmouse <- rownames(inputdatgseamousedown)
