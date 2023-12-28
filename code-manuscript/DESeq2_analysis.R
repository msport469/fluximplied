## ----include=FALSE---------------------------------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(tibble)

library(DESeq2)
library(ashr)

library(DOSE)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(ggplot2)
library(viridis)
library(enrichR)
library(gridExtra)


## --------------------------------------------------------------------------------------------------------------------------


setwd("/Users/adam/Box/GeberA_TophamLab/03-GTEx_Human_Metabolism/Analyses")
attr <- read.table("~/Box/GeberA_TophamLab/03-GTEx_Human_Metabolism/metadata/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header = TRUE, row.names = 1, sep = "\t", quote = "", stringsAsFactors = FALSE, fill = FALSE) # this file was downloaded from version 8 of bulk RNAseq from GTEx. If you can't find it please contact us.

attr_filter <- filter(attr, SMAFRZE == 'RNASEQ')
attr_filter$SUBJID <- unlist(lapply(strsplit(rownames(attr_filter), '-'), FUN=function(x){paste(x[1],x[2],sep="-")}))

attr_adipose <- filter(attr_filter, SMTSD == 'Adipose - Subcutaneous')
attr_liver <- filter(attr_filter, SMTSD == 'Liver')
attr_putamen <- filter(attr_filter, SMTSD == 'Brain - Putamen (basal ganglia)')

pheno <- read.table("~/Box/GeberA_TophamLab/03-GTEx_Human_Metabolism/metadata/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep = "\t", quote = "", fill = TRUE, header = TRUE) # this file was downloaded from version 8 of bulk RNAseq from GTEx. If you can't find it please contact us.

pheno_filter <- filter(pheno, DTHHRDY == 1) # restricting to a Hardy Scale value of 1 means that the process of dying should, in theory, have a more minimal effect on altering tissue metabolism

adipose <-subset(attr_adipose, attr_adipose$SUBJID %in% pheno_filter$SUBJID) # retains only those SUBJIDs that match pheno_filter
liver <-subset(attr_liver, attr_liver$SUBJID %in% pheno_filter$SUBJID)
putamen <-subset(attr_putamen, attr_putamen$SUBJID %in% pheno_filter$SUBJID)


## --------------------------------------------------------------------------------------------------------------------------
counts <- fread("~/Box/GeberA_TophamLab/03-GTEx_Human_Metabolism/countData/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", verbose = TRUE) # works with 247s total runtime # this file was downloaded from version 8 of bulk RNAseq from GTEx. It's probably easier for you to use the "tissueCounts.csv file we generated after subsetting on relevant data (see a few lines below). If you can't find it please contact us.


## --------------------------------------------------------------------------------------------------------------------------
adiposeCounts <- subset(counts, select = c("Name", rownames(adipose)))
liverCounts <-  subset(counts, select = c("Name", rownames(liver)))
putamenCounts <- subset(counts, select = c("Name", rownames(putamen)))
tissueCounts <- bind_cols(adiposeCounts, liverCounts[,-(1)], putamenCounts[,-(1)]) # exclude the "Name" column from the second two dataframes in order to avoid duplication
write.csv(tissueCounts, file = "tissueCounts.csv")


## --------------------------------------------------------------------------------------------------------------------------
metadata <- bind_rows(adipose, liver, putamen)
metadata <- metadata %>%
  dplyr::select(SUBJID, SMTSD)


## --------------------------------------------------------------------------------------------------------------------------
tissueCounts <- column_to_rownames(tissueCounts, var = "Name")
all(rownames(metadata) == colnames(tissueCounts)) # looks like we have matching between columns and rows, as DESeq2 requires
dds <- DESeqDataSetFromMatrix(countData = tissueCounts, colData = metadata, design = ~SUBJID + SMTSD)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds1 <- nbinomWaldTest(dds, maxit=500) # per Michael Love, this may be necessary to address the 15 rows that didn't converge. bumping up the max iterations to 500 still leaves 7 rows that don't converge, though.
ddsClean <- dds1[which(mcols(dds1)$betaConv),]
res <- results(ddsClean)

## --------------------------------------------------------------------------------------------------------------------------
library(gplots, RColorBrewer)
vsd <- vst(ddsClean)
head(assay(vsd))

par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(ddsClean, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(vsd)[, 1:2], col="#00000020", pch=20, cex=0.3 )

sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$SMTSD)
colnames(sampleDistMatrix) <- paste( vsd$SMTSD)
heatmap.2( sampleDistMatrix, trace="none", margins = c(8,10)) # default palette

plotPCA(vsd, intgroup = "SMTSD")


## --------------------------------------------------------------------------------------------------------------------------
AdiposevLiver.lfc <- lfcShrink(ddsClean, contrast = c("SMTSD", "Adipose - Subcutaneous", "Liver"), type = "ashr")
AdiposevPutamen.lfc <- lfcShrink(ddsClean, contrast = c("SMTSD", "Adipose - Subcutaneous", "Brain - Putamen (basal ganglia)"), type = "ashr")
LivervPutamen.lfc <- lfcShrink(ddsClean, contrast = c("SMTSD", "Liver", "Brain - Putamen (basal ganglia)"), type = "ashr")

AvL.DE <- as.data.frame(AdiposevLiver.lfc[which(AdiposevLiver.lfc$padj<0.05 & AdiposevLiver.lfc$baseMean > 10 & abs(AdiposevLiver.lfc$log2FoldChange)>0.5),])
AvP.DE <- as.data.frame(AdiposevPutamen.lfc[which(AdiposevPutamen.lfc$padj<0.05 & AdiposevLiver.lfc$baseMean > 10 & abs(AdiposevPutamen.lfc$log2FoldChange)>0.5),])
LvP.DE <- as.data.frame(LivervPutamen.lfc[which(LivervPutamen.lfc$padj<0.05 & AdiposevLiver.lfc$baseMean > 10 & abs(LivervPutamen.lfc$log2FoldChange)>0.5),])

rownames(AvL.DE) <- do.call(rbind, strsplit(as.character(rownames(AvL.DE)), "\\."))[,1] # remove versioning suffix from each Ensembl ID
rownames(AvP.DE) <- do.call(rbind, strsplit(as.character(rownames(AvP.DE)), "\\."))[,1]
rownames(LvP.DE) <- do.call(rbind, strsplit(as.character(rownames(LvP.DE)), "\\."))[,1]

AvL.DE$SYMBOL <- mapIds(org.Hs.eg.db, keys = rownames(AvL.DE), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "asNA")
AvP.DE$SYMBOL <- mapIds(org.Hs.eg.db, keys = rownames(AvP.DE), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "asNA")
LvP.DE$SYMBOL <- mapIds(org.Hs.eg.db, keys = rownames(LvP.DE), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "asNA")

AvL.DE <- na.omit(AvL.DE)
AvP.DE <- na.omit(AvP.DE)
LvP.DE <- na.omit(LvP.DE)

AvL.up <- AvL.DE[which(AvL.DE$log2FoldChange > 0),] # separate by LFC sign (i.e. per sample) for each pairwise comparison
AvL.down <- AvL.DE[which(AvL.DE$log2FoldChange < 0),]
AvP.up <- AvP.DE[which(AvP.DE$log2FoldChange > 0),]
AvP.down <- AvP.DE[which(AvP.DE$log2FoldChange < 0),]
LvP.up <- LvP.DE[which(LvP.DE$log2FoldChange > 0),]
LvP.down <- LvP.DE[which(LvP.DE$log2FoldChange < 0),]

rownames(AvL.up) <- AvL.up$SYMBOL # all unique values
AvL.down <- AvL.down[!rownames(AvL.down) %in% c("ENSG00000271858", "ENSG00000261186"),] # retain only rows that don't match these specific rownames, as detailed below
rownames(AvL.down) <- AvL.down$SYMBOL # based on non-unique values, I manually removed the rows corresponding to non-coding transcripts mapping to gene symbols 'CYB561D2', 'LINC01238'
rownames(AvP.up) <- AvP.up$SYMBOL # all unique values
AvP.down <- AvP.down[!rownames(AvP.down) %in% c("ENSG00000243902", "ENSG00000204334"),]
rownames(AvP.down) <- AvP.down$SYMBOL # based on non-unique values, I manually removed the rows corresponding to non-coding transcripts mapping to gene symbols ‘ELFN2’, ‘ERICH2’
LvP.up <- LvP.up[!rownames(LvP.up) %in% ("ENSG00000174353"),]
rownames(LvP.up) <- LvP.up$SYMBOL # based on non-unique values, I manually removed the rows corresponding to pseudogenes mapping to gene symbols ‘TRIM74’
LvP.down <- LvP.down[!rownames(LvP.down) %in% ("ENSG00000243902"),]
rownames(LvP.down) <- LvP.down$SYMBOL # based on non-unique values, I manually removed the rows corresponding to non-coding transcripts mapping to gene symbols ‘ELFN2’

AvL.combined <- rbind(AvL.down, AvL.up) # recombine the up/down datasets for input into fluximplied
AvP.combined <- rbind(AvP.down, AvP.up)
LvP.combined <- rbind(LvP.down, LvP.up)

setwd("~/Documents/GitHub/fluximplieddev")

saveRDS(AvL.combined, "humandata/AdiposevLiver.RDS")
saveRDS(AvP.combined, "humandata/AdiposevPutamen.RDS")
saveRDS(LvP.combined, "humandata/LivervPutamen.RDS")


## --------------------------------------------------------------------------------------------------------------------------
fluximplied(AvL.combined, species = "Hsa", geneform = "symbol", inputformat = "df", padjcolname = "padj", pcutoff = 0.05)
fluximplied(AvP.combined, species = "Hsa", geneform = "symbol", inputformat = "df", padjcolname = "padj", pcutoff = 0.05)
fluximplied(LvP.combined, species = "Hsa", geneform = "symbol", inputformat = "df", padjcolname = "padj", pcutoff = 0.05)

#length(rownames(AvL.DE)) == length(unique(rownames(AvL.DE))) # okay, so this lets us know that with the GTEx ENSEMBL ID format (which includes versioning info after the decimal) is at least unique per row

AvL.up.results <- enrichr(AvL.up, databases = "Reactome_2016")
AvL.down.results <- enrichr(AvL.down, databases = "Reactome_2016")
plotEnrich(AvL.up.results[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(AvL.down.results[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

AvP.up.results <- enrichr(AvP.up, databases = "Reactome_2016")
AvP.down.results <- enrichr(AvP.down, databases = "Reactome_2016")
plotEnrich(AvP.up.results[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(AvP.down.results[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

LvP.up.results <- enrichr(LvP.up, databases = "Reactome_2016")
LvP.down.results <- enrichr(LvP.down, databases = "Reactome_2016")
plotEnrich(LvP.up.results[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
plotEnrich(LvP.down.results[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

sessionInfo()
#citations, using a method as described in https://stackoverflow.com/questions/27535628/how-do-i-tell-which-r-packages-to-cite-in-my-paper 
packages_in_use <- c(sessionInfo()$basePkgs, names( sessionInfo()$loadedOnly ) )
the_citations_list <- lapply( X=packages_in_use, FUN=citation)
sink(file='the_citation_list.txt')
the_citations_list
sink(file=NULL)
####### FIN #######
