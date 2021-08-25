#load libs
library(ggplot2)
library(viridis)
library(enrichR)
library(gridExtra)
library(tidyr)
library(fluximplied)
library(cowplot)
#read in humans
AdiposevLiver<-readRDS('humandata/AdiposevLiver.RDS')
AdiposevPutamen<-readRDS('humandata/AdiposevPutamen.RDS')
LivervPutamen<-readRDS('humandata/LivervPutamen.RDS')

#LFC cutoff
inputdatgseaAvLup <- inputdatgseaAvL[inputdatgseaAvL[,"log2FoldChange"]>0.5,]
inputdatgseaAvLdown <- inputdatgseaAvL[inputdatgseaAvL[,"log2FoldChange"]<(-0.5),]
upgenesAvL<-rownames(inputdatgseaAvLup)
downgenesAvL<-rownames(inputdatgseaAvLdown)

inputdatgseaAvPup <- inputdatgseaAvP[inputdatgseaAvP[,"log2FoldChange"]>0.5,]
inputdatgseaAvPdown <- inputdatgseaAvP[inputdatgseaAvP[,"log2FoldChange"]<(-0.5),]
upgenesAvP<-rownames(inputdatgseaAvPup)
downgenesAvP<-rownames(inputdatgseaAvPdown)

inputdatgseaLvPup <- inputdatgseaLvP[inputdatgseaLvP[,"log2FoldChange"]>0.5,]
inputdatgseaLvPdown <- inputdatgseaLvP[inputdatgseaLvP[,"log2FoldChange"]<(-0.5),]
upgenesLvP<-rownames(inputdatgseaLvPup)
downgenesLvP<-rownames(inputdatgseaLvPdown)

####read in MOUSEYMOUSEYMOUSEY YEAH####
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

#set up dbs for gsea
dbs <- c("KEGG_2019_Human",'Reactome_2016')

##### Mouse #####
##### Mouse #####
##### Mouse #####
#get ready to plot 
list_up <- c(upgenes)
list_down <- c(downgenes)
eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)



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
gos$Term <- sub("\ Homo.*$", "", gos$Term) # remove Reactome-specific IDs
gos$Term <- sub("independent", "\nindependent", gos$Term) # insert linebreak
reactmouse<-ggplot(gos, aes(x=reorder(Term,Overlap), y=Overlap , label=Overlap)) + 
  geom_bar(stat='identity', aes(fill=logpadj), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Overlap') +
  ylim(-1, 1) +
  xlab(NULL)+
  labs(fill = bquote(-Log(P[adj]))) +
  coord_flip()
reactmouse

write.csv(eup$Reactome_2016,'enrichrcsvs/upreactome-mouse-LFC_0.5_.csv')
write.csv(edown$Reactome_2016,'enrichrcsvs/downreactome-mouse-LFC_0.5_.csv')
###### Hoooman #####
###############
#############
###### Hoooman #####

#AvL
list_up <- c(upgenesAvL)
list_down <- c(downgenesAvL)
eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)

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
AvLkeggplot<-ggplot(gos, aes(x=reorder(Term,Overlap), y=Overlap , label=Overlap)) + 
  geom_bar(stat='identity', aes(fill=logpadj), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Overlap') +
  ylim(-1, 1) +
  xlab(NULL)+
  labs(fill = bquote(-Log(P[adj]))) +
  coord_flip()
AvLkeggplot

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
gos$Term <- sub("\ Homo.*$", "", gos$Term) # remove Reactome-specific IDs
gos$Term <- sub("cell\ spreading", "\ncell\ spreading", gos$Term) # insert linebreak
AvLreactplot<-ggplot(gos, aes(x=reorder(Term,Overlap), y=Overlap , label=Overlap)) + 
  geom_bar(stat='identity', aes(fill=logpadj), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Overlap') +
  ylim(-1, 1) +
  xlab(NULL)+
  labs(fill = bquote(-Log(P[adj]))) +
  coord_flip()
AvLreactplot

##### write human csvs #####
write.csv(eup$KEGG_2019_Human,'./enrichrcsvs/AvL-upKEGG_2019_Human-LFC_0.5_.csv')
write.csv(eup$Reactome_2016,'./enrichrcsvs/AvL-upReactome_2016_Human-LFC_0.5_.csv')

write.csv(edown$KEGG_2019_Human,'./enrichrcsvs/AvL-downKEGG_2019_Human-LFC_0.5_.csv')
write.csv(edown$Reactome_2016,'./enrichrcsvs/AvL-downReactome_2016_Human-LFC_0.5_.csv')


#AvP
list_up <- c(upgenesAvP)
list_down <- c(downgenesAvP)
eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)

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
AvPkeggplot<-ggplot(gos, aes(x=reorder(Term,Overlap), y=Overlap , label=Overlap)) + 
  geom_bar(stat='identity', aes(fill=logpadj), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Overlap') +
  ylim(-1, 1) +
  xlab(NULL)+
  labs(fill = bquote(-Log(P[adj]))) +
  coord_flip()
AvPkeggplot

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
gos$Term <- sub("\ Homo.*$", "", gos$Term) # remove Reactome-specific IDs
gos$Term <- sub("Triggers", "\nTriggers", gos$Term) # insert linebreak
AvPreactplot<-ggplot(gos, aes(x=reorder(Term,Overlap), y=Overlap , label=Overlap)) + 
  geom_bar(stat='identity', aes(fill=logpadj), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Overlap') +
  ylim(-1, 1) +
  xlab(NULL)+
  labs(fill = bquote(-Log(P[adj]))) +
  coord_flip()
AvPreactplot

##### write human csvs #####
write.csv(eup$KEGG_2019_Human,'./enrichrcsvs/AvP-upKEGG_2019_Human-LFC_0.5_.csv')
write.csv(eup$Reactome_2016,'./enrichrcsvs/AvP-upReactome_2016_Human-LFC_0.5_.csv')

write.csv(edown$KEGG_2019_Human,'./enrichrcsvs/AvP-downKEGG_2019_Human-LFC_0.5_.csv')
write.csv(edown$Reactome_2016,'./enrichrcsvs/AvP-downReactome_2016_Human-LFC_0.5_.csv')

#LvP
list_up <- c(upgenesLvP)
list_down <- c(downgenesLvP)
eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)

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
LvPkeggplot<-ggplot(gos, aes(x=reorder(Term,Overlap), y=Overlap , label=Overlap)) + 
  geom_bar(stat='identity', aes(fill=logpadj), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Overlap') +
  ylim(-1, 1) +
  xlab(NULL)+
  labs(fill = bquote(-Log(P[adj]))) +
  coord_flip()
LvPkeggplot

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
gos$Term <- sub("\ Homo.*$", "", gos$Term) # remove Reactome-specific IDs
gos$Term <- sub("Triggers", "\nTriggers", gos$Term) # insert linebreak
LvPreactplot<-ggplot(gos, aes(x=reorder(Term,Overlap), y=Overlap , label=Overlap)) + 
  geom_bar(stat='identity', aes(fill=logpadj), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  # labs(title= "") +
  ylab('Overlap') +
  ylim(-1, 1) +
  xlab(NULL)+
  labs(fill = bquote(-Log(P[adj]))) +
  coord_flip()
LvPreactplot

##### write human csvs #####
write.csv(eup$KEGG_2019_Human,'./enrichrcsvs/LvP-upKEGG_2019_Human-LFC_0.5_.csv')
write.csv(eup$Reactome_2016,'./enrichrcsvs/LvP-upReactome_2016_Human-LFC_0.5_.csv')

write.csv(edown$KEGG_2019_Human,'./enrichrcsvs/LvP-downKEGG_2019_Human-LFC_0.5_.csv')
write.csv(edown$Reactome_2016,'./enrichrcsvs/LvP-downReactome_2016_Human-LFC_0.5_.csv')

##### fluximplied stuff #####
##### fluximplied stuff #####
##### fluximplied stuff #####
fluximplied(inputdatgseamouse,padjcolname = 'weighted_pvalue')
mousefluxplot<-fluximpliedplot

fluximplied(AdiposevLiver,species = 'hsa',padjcolname = 'padj')
AdiposevLiverfluxplot<-fluximpliedplot

fluximplied(AdiposevPutamen,species = 'hsa',padjcolname = 'padj')
AdiposevPutamenfluxplot<-fluximpliedplot

fluximplied(LivervPutamen,species = 'hsa',padjcolname = 'padj')
LivervPutamenfluxplot<-fluximpliedplot

#### gridarranging plots ####
#kegg
a<-AvLkeggplot+ggtitle('Adipose vs. Liver')
s<-AvPkeggplot+ggtitle('Adipose vs. Putamen')
d<-LvPkeggplot+ggtitle('Liver vs. Putamen')
plot_grid(a,s,d, ncol = 1, align = "v")

grid.arrange(a,s,d,ncol=1)

#react
z<-AvLreactplot+ggtitle('Adipose vs. Liver')
x<-AvPreactplot+ggtitle('Adipose vs. Putamen')
c<-LvPreactplot+ggtitle('Liver vs. Putamen')
v<-reactmouse+ggtitle('TRM vs circulating T cells')
plot_grid(z,x,c,v, ncol = 1, align = "v")

grid.arrange(z,x,c,v,ncol=1)

#fluximplied
q<-AdiposevLiverfluxplot+ggtitle('Adipose vs. Liver')
w<-AdiposevPutamenfluxplot+ggtitle('Adipose vs. Putamen')
e<-LivervPutamenfluxplot+ggtitle('Liver vs. Putamen')
r<-mousefluxplot+ggtitle('TRM vs circulating T cells')
grid.arrange(q,w,e,r,ncol=1)
plot_grid(q,w,e,r, ncol = 1, align = "v")

# Adipose v Liver
plot_grid(a, z, q, ncol = 1, align = "v")

######
grid.arrange(a,z,q,ncol=1)
grid.arrange(arrangeGrob(a,z),q,nrow=3)
#AvLplots<-grid.arrange(arrangeGrob(a,z,q,nrow=2,widths=c(0.5,0.5,1)))


#ggsave('humanbarplots.png',plot=barplots,dpi=600,path='./compare/',units = 'cm',width=16.6,height=6,scale=2)




