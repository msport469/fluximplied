library(ggplot2)
library(viridis)
library(enrichR)

#load in deseqresult
exampledeseqresultdataframe<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/exampledeseqresultdataframe.csv',row.names = c(1))
inputdat=exampledeseqresultdataframe
inputdatgsea<-inputdat
inputdatgsea <- inputdatgsea[inputdatgsea[,"weighted_pvalue"]<0.05,]
inputdatgseaup <- inputdatgsea[inputdatgsea[,"log2FoldChange"]>0,]
inputdatgseadown <- inputdatgsea[inputdatgsea[,"log2FoldChange"]<0,]


upgenes <- rownames(inputdatgseaup)
downgenes <- rownames(inputdatgseadown)


#Enrichr
listEnrichrDbs()
alldbs <- listEnrichrDbs()
dbs <- c("KEGG_2019_Mouse","KEGG_2019_Human",'Reactome_2016')
theme_set(theme_bw())
list_up <- c(upgenes)
list_down <- c(downgenes)
eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)

#KEGG Mouse
up <- eup$KEGG_2019_Mouse
down <- edown$KEGG_2019_Mouse
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="Combined scores", 
       title= "GSEA KEGG_2019_Mouse") + 
  coord_flip()

#KEGG Human
up <- eup$KEGG_2019_Human
down <- edown$KEGG_2019_Human
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="Combined scores", 
       title= "GSEA KEGG_2019_Human") + 
  coord_flip()


#Reactome 2016
up <- eup$Reactome_2016
down <- edown$Reactome_2016
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="Combined scores", 
       title= "GSEA Reactome_2016") + 
  coord_flip()

write.csv(eup$KEGG_2019_Mouse,'./compare/upKEGG_2019_Mouse.csv')
write.csv(eup$KEGG_2019_Human,'./compare/upKEGG_2019_Human.csv')
write.csv(eup$Reactome_2016,'./compare/upReactome_2016.csv')

write.csv(edown$KEGG_2019_Mouse,'./compare/downKEGG_2019_Mouse.csv')
write.csv(edown$KEGG_2019_Human,'./compare/downKEGG_2019_Human.csv')
write.csv(edown$Reactome_2016,'./compare/downReactome_2016.csv')


#now make the same things but with lfc more stringent
inputdatgseaup <- inputdatgsea[inputdatgsea[,"log2FoldChange"]>1,]
inputdatgseadown <- inputdatgsea[inputdatgsea[,"log2FoldChange"]<(-1),]
upgenes <- rownames(inputdatgseaup)
downgenes <- rownames(inputdatgseadown)


#Enrichr
listEnrichrDbs()
alldbs <- listEnrichrDbs()
dbs <- c("KEGG_2019_Mouse","KEGG_2019_Human",'Reactome_2016')
theme_set(theme_bw())
list_up <- c(upgenes)
list_down <- c(downgenes)
eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)

#KEGG Mouse
up <- eup$KEGG_2019_Mouse
down <- edown$KEGG_2019_Mouse
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="Combined scores", 
       title= "GSEA KEGG_2019_Mouse") + 
  coord_flip()

#KEGG Human
up <- eup$KEGG_2019_Human
down <- edown$KEGG_2019_Human
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="Combined scores", 
       title= "GSEA KEGG_2019_Human") + 
  coord_flip()


#Reactome 2016
up <- eup$Reactome_2016
down <- edown$Reactome_2016
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="Combined scores", 
       title= "GSEA Reactome_2016") + 
  coord_flip()

write.csv(eup$KEGG_2019_Mouse,'./compare/upKEGG_2019_Mouse-stringentlfc.csv')
write.csv(eup$KEGG_2019_Human,'./compare/upKEGG_2019_Human-stringentlfc.csv')
write.csv(eup$Reactome_2016,'./compare/upReactome_2016-stringentlfc.csv')

write.csv(edown$KEGG_2019_Mouse,'./compare/downKEGG_2019_Mouse-stringentlfc.csv')
write.csv(edown$KEGG_2019_Human,'./compare/downKEGG_2019_Human-stringentlfc.csv')
write.csv(edown$Reactome_2016,'./compare/downReactome_2016-stringentlfc.csv')
