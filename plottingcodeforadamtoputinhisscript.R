library(ggplot2)
library(viridis)
library(enrichR)
library(gridExtra)
library(tidyr)
#Enrichr
listEnrichrDbs()
alldbs <- listEnrichrDbs()
dbs <- c("KEGG_2019_Human",'Reactome_2016')
theme_set(theme_grey())
list_up <- c(rownames(AvL.up))
list_down <- c(rownames(AvL.down))
eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)

#KEGG Mouse
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

#Reactome 2016
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



write.csv(eup$KEGG_2019_Human,'./compare/upKEGG_2019_Human-LFC_0.5_.csv')
write.csv(eup$Reactome_2016,'./compare/upReactome_2016-LFC_0.5_.csv')

write.csv(edown$KEGG_2019_Human,'./compare/downKEGG_2019_Human-LFC_0.5_.csv')
write.csv(edown$Reactome_2016,'./compare/downReactome_2016-LFC_0.5_.csv')

barplots<-grid.arrange(arrangeGrob(reactomebarplot,keggbarplot,nrow=1,widths=c(1.5,1)))
barplots
ggsave('humanbarplots.png',plot=barplots,dpi=600,path='./compare/',units = 'cm',width=16.6,height=6,scale=2)