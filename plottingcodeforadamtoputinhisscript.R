#load libs
library(ggplot2)
library(viridis)
library(enrichR)
library(gridExtra)
library(tidyr)

#Enrichr
#set up dbs
dbs <- c("KEGG_2019_Human",'Reactome_2016')
#read in humans
AdiposevLiver<-readRDS('humandata/AdiposevLiver.RDS')
AdiposevPutamen<-readRDS('humandata/AdiposevPutamen.RDS')
LivervPutamen<-readRDS('humandata/LivervPutamen.RDS')
#p value subset
inputdatgseaAvL <- AdiposevLiver[AdiposevLiver[,"weighted_pvalue"]<0.05,]
inputdatgseaAvP <- AdiposevPutamen[AdiposevPutamen[,"weighted_pvalue"]<0.05,]
inputdatgseaLvP <- LivervPutamen[LivervPutamen[,"weighted_pvalue"]<0.05,]
#LFC cutoff
inputdatgseaAvLup <- inputdatgseaAvL[inputdatgseaAvL[,"log2FoldChange"]>0.5,]
inputdatgseaAvLdown <- inputdatgseaAvL[inputdatgseaAvL[,"log2FoldChange"]<(-0.5),]

inputdatgseaAvPup <- inputdatgseaAvP[inputdatgseaAvP[,"log2FoldChange"]>0.5,]
inputdatgseaAvPdown <- inputdatgseaAvP[inputdatgseaAvP[,"log2FoldChange"]<(-0.5),]

inputdatgseaLvPup <- inputdatgseaLvP[inputdatgseaLvP[,"log2FoldChange"]>0.5,]
inputdatgseaLvPdown <- inputdatgseaLvP[inputdatgseaLvP[,"log2FoldChange"]<(-0.5),]

####MOUSEYMOUSEYMOUSEY YEAH####
inputdatgseamouse<-readRDS('mousedata/DPDNdfmouse.RDS')
inputdatgseamouse$gene<-rownames(inputdatgseamouse)
#subset on p values
inputdatgseamouse <- inputdatgseamouse[inputdatgseamouse[,"weighted_pvalue"]<0.05,]
#differentiate up and down
inputdatgseamouseup <- inputdatgseamouse[inputdatgseamouse[,"log2FoldChange"]>0.5,]
inputdatgseamousedown <- inputdatgseamouse[inputdatgseamouse[,"log2FoldChange"]<(-0.5),]
#save genes as vectors.
upgenesmouse <- rownames(inputdatgseamouseup)
downgenesmouse <- rownames(inputdatgseamousedown)
upgenes<-upgenesmouse
downgenes<-downgenesmouse

#get ready to plot 
list_up <- c(upgenes)
list_down <- c(downgenes)
eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)

#KEGG Mouse
#####there are no up or down pathways in KEGG mouse cutoff lfc=0.5 or -0.5
# up <- eup$KEGG_2019_Mouse
# down <- edown$KEGG_2019_Mouse
# up$type <- "up"
# down$type <- "down"
# up <- up[up$Adjusted.P.value<.05,]
# up <- up[order(up$Combined.Score), ]  # sort
# down <- down[down$Adjusted.P.value<0.05,]
# down <- down[order(down$Combined.Score), ]  # sort
# down$Combined.Score <- (-1) * down$Combined.Score
# gos <- rbind(down,up)
# gos <- na.omit(gos) # Diverging Barcharts
# keggbarplot<-ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
#   geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
#   scale_fill_viridis() + 
#   # labs(title= "") +
#   ylab('Combined Score') +
#   xlab(NULL)+
#   labs(fill = bquote('P'['adj']))+
#   coord_flip()


#Reactome 2016 mouse
up <- eup$Reactome_2016
down <- edown$Reactome_2016
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score, decreasing = TRUE), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score, decreasing = TRUE), ]  # sort
up<-up[1:10,]
up <- separate(up, Overlap, sep = "/", into = c("Num", "Denom"))
up$Num <- as.numeric(up$Num)
up$Denom <- as.numeric(up$Denom)
up$Overlap <- up$Num / up$Denom
down<-down[1:10,]
down <- separate(down, Overlap, sep = "/", into = c("Num", "Denom"))
down$Num <- as.numeric(down$Num)
down$Denom <- as.numeric(down$Denom)
down$Overlap <- down$Num / down$Denom
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
gos$logpadj <- -log10(gos$Adjusted.P.value)
gos$Overlap[1:10] <- -(gos$Overlap[1:10])
ggplot(gos, aes(x=reorder(Term,Overlap), y=Overlap , label=Overlap)) + 
  geom_bar(stat='identity', aes(fill=logpadj), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Overlap') +
  ylim(-1, 1) +
  xlab(NULL)+
  labs(fill = bquote('-log P'['adj']))+
  coord_flip()
reactomebarplot








#we need to add lines of code to read in files. 
#adam, can you please put .RDS files in subfolder of fluximplied dev called
#human data
readRDS('humandata/AdiposevLiver.RDS')

list_up <- c(upgenes)
list_down <- c(downgenes)

eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)

#KEGG Hoooman
up <- eup$KEGG_2019_Human
down <- edown$KEGG_2019_Human
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score, decreasing = TRUE), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score, decreasing = TRUE), ]  # sort
up<-up[1:10,]
up <- separate(up, Overlap, sep = "/", into = c("Num", "Denom"))
up$Num <- as.numeric(up$Num)
up$Denom <- as.numeric(up$Denom)
up$Overlap <- up$Num / up$Denom
down<-down[1:10,]
down <- separate(down, Overlap, sep = "/", into = c("Num", "Denom"))
down$Num <- as.numeric(down$Num)
down$Denom <- as.numeric(down$Denom)
down$Overlap <- down$Num / down$Denom
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
gos$logpadj <- -log10(gos$Adjusted.P.value)
gos$Overlap[1:10] <- -(gos$Overlap[1:10])
ggplot(gos, aes(x=reorder(Term,Overlap), y=Overlap , label=Overlap)) + 
  geom_bar(stat='identity', aes(fill=logpadj), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Overlap') +
  ylim(-1, 1) +
  xlab(NULL)+
  labs(fill = bquote('-log P'['adj']))+
  coord_flip()

temp <- up
temp <- separate(temp, Overlap, sep = "/", into = c("Num", "Denom"))
temp1 <- temp
temp1$Num <- as.numeric(temp$Num)
temp1$Denom <- as.numeric(temp$Denom)
temp1$Overlap <- temp1$Num / temp1$Denom

#Reactome 2016 human
up <- eup$Reactome_2016
down <- edown$Reactome_2016
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score, decreasing = TRUE), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score, decreasing = TRUE), ]  # sort
up<-up[1:10,]
up <- separate(up, Overlap, sep = "/", into = c("Num", "Denom"))
up$Num <- as.numeric(up$Num)
up$Denom <- as.numeric(up$Denom)
up$Overlap <- up$Num / up$Denom
down<-down[1:10,]
down <- separate(down, Overlap, sep = "/", into = c("Num", "Denom"))
down$Num <- as.numeric(down$Num)
down$Denom <- as.numeric(down$Denom)
down$Overlap <- down$Num / down$Denom
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
gos$logpadj <- -log10(gos$Adjusted.P.value)
gos$Overlap[1:10] <- -(gos$Overlap[1:10])
ggplot(gos, aes(x=reorder(Term,Overlap), y=Overlap , label=Overlap)) + 
  geom_bar(stat='identity', aes(fill=logpadj), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Overlap') +
  ylim(-1, 1) +
  xlab(NULL)+
  labs(fill = bquote('-log P'['adj']))+
  coord_flip()
reactomebarplot











write.csv(eup$KEGG_2019_Human,'./compare/upKEGG_2019_Human-LFC_0.5_.csv')
write.csv(eup$Reactome_2016,'./compare/upReactome_2016-LFC_0.5_.csv')

write.csv(edown$KEGG_2019_Human,'./compare/downKEGG_2019_Human-LFC_0.5_.csv')
write.csv(edown$Reactome_2016,'./compare/downReactome_2016-LFC_0.5_.csv')

barplots<-grid.arrange(arrangeGrob(reactomebarplot,keggbarplot,nrow=1,widths=c(1.5,1)))
barplots
ggsave('humanbarplots.png',plot=barplots,dpi=600,path='./compare/',units = 'cm',width=16.6,height=6,scale=2)