library(ggplot2)
library(viridis)
library(enrichR)
library(gridExtra)
#load in deseqresult
exampledeseqresultdataframe<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/exampledeseqresultdataframe.csv',row.names = c(1))
inputdat=exampledeseqresultdataframe
inputdatgsea<-inputdat
inputdat$gene<-rownames(inputdat)
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
inputdatgseaup <- inputdatgsea[inputdatgsea[,"log2FoldChange"]>0.5,]
inputdatgseadown <- inputdatgsea[inputdatgsea[,"log2FoldChange"]<(-0.5),]
upgenes <- rownames(inputdatgseaup)
downgenes <- rownames(inputdatgseadown)


#Enrichr
listEnrichrDbs()
alldbs <- listEnrichrDbs()
dbs <- c("KEGG_2019_Mouse","KEGG_2019_Human",'Reactome_2016')
#theme_set(theme_bw())
theme_set(theme_grey())
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
keggbarplot<-ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Combined Score') +
  xlab(NULL)+
  labs(fill = bquote('P'['adj']))+
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
up<-up[1:5,]
down<-down[1:5,]
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
reactomebarplot<-ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
 # labs(title= "") +
  ylab('Combined Score') +
  xlab(NULL)+
  labs(fill = bquote('P'['adj']))+
  coord_flip()


write.csv(eup$KEGG_2019_Mouse,'./compare/upKEGG_2019_Mouse-LFC_0.5_.csv')
write.csv(eup$KEGG_2019_Human,'./compare/upKEGG_2019_Human-LFC_0.5_.csv')
write.csv(eup$Reactome_2016,'./compare/upReactome_2016-LFC_0.5_.csv')

write.csv(edown$KEGG_2019_Mouse,'./compare/downKEGG_2019_Mouse-LFC_0.5_.csv')
write.csv(edown$KEGG_2019_Human,'./compare/downKEGG_2019_Human-LFC_0.5_.csv')
write.csv(edown$Reactome_2016,'./compare/downReactome_2016-LFC_0.5_.csv')

barplots<-grid.arrange(arrangeGrob(reactomebarplot,keggbarplot,nrow=1,widths=c(1.5,1)))
barplots
ggsave('mousebarplots.png',plot=barplots,dpi=600,path='./compare/',units = 'cm',width=16.6,height=6,scale=2)
