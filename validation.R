dp49a<-read.csv('deseq2result.dpvscd49a.csv') 
head(dp49a)
dim(dp49a)
updp49a<-subset(dp49a,dp49a$log2FoldChange>0)
head(updp49a)
dim(updp49a)
downdp49a<-subset(dp49a,dp49a$log2FoldChange<0)
head(downdp49a)
dim(downdp49a)

names(dp49a)[1] <- "genes"
head(dp49a)
sub<-subset(dp49a, dp49a$genes %in% RLSdatabase$mouse.gene.symbol)
sub



input<-read.csv('CD103vsDN.all.csv',row.names = c(1))
input$gene<-rownames(input)
input1=subset(input,input$log2FoldChange>1)#
input2=subset(input,input$log2FoldChange<(-1))#

input<-read.csv('CD103vsCD49a.all.csv',row.names = c(1))
input$gene<-rownames(input)
input3=subset(input,input$log2FoldChange>1)#
input4=subset(input,input$log2FoldChange<(-1))#

input<-read.csv('CD103vsDP.all.csv',row.names = c(1))
input$gene<-rownames(input)
input5=subset(input,input$log2FoldChange>1)#
input6=subset(input,input$log2FoldChange<(-1))

input<-read.csv('CD49avsDN.all.csv',row.names = c(1))
input$gene<-rownames(input)
input7=subset(input,input$log2FoldChange>1)
input8=subset(input,input$log2FoldChange<(-1))

input<-read.csv('CD49avsDP.all.csv',row.names = c(1))
input$gene<-rownames(input)
input9=subset(input,input$log2FoldChange>1)
input10=subset(input,input$log2FoldChange<(-1))

input<-read.csv('DPvsDN.all.csv',row.names = c(1))
input$gene<-rownames(input)
input11=subset(input,input$log2FoldChange>1)
input12=subset(input,input$log2FoldChange<(-1))

inputlist<-list(inputlist)
inputlist<-paste0('input',1:12)




flux1<-fluximplied(input1,
                   species,
                   geneformat,
                   inputformat,
                   padjcolname)
flux2<-fluximplied(input2,
                   species,
                   geneformat,
                   inputformat,
                   padjcolname)
flux3<-fluximplied(input3,
                   species,
                   geneformat,
                   inputformat,
                   padjcolname)
flux4<-fluximplied(input4,
                   species,
                   geneformat,
                   inputformat,
                   padjcolname)
flux5<-fluximplied(input5,
                   species,
                   geneformat,
                   inputformat,
                   padjcolname)
flux6<-fluximplied(input6,
                   species,
                   geneformat,
                   inputformat,
                   padjcolname)
flux7<-fluximplied(input7,
                   species,
                   geneformat,
                   inputformat,
                   padjcolname)
flux8<-fluximplied(input8,
                   species,
                   geneformat,
                   inputformat,
                   padjcolname)
flux9<-fluximplied(input9,
                   species,
                   geneformat,
                   inputformat,
                   padjcolname)
flux10<-fluximplied(input10,
                    species,
                    geneformat,
                    inputformat,
                    padjcolname)
flux11<-fluximplied(input11,
                    species,
                    geneformat,
                    inputformat,
                    padjcolname)
flux12<-fluximplied(input12,
                    species,
                    geneformat,
                    inputformat,
                    padjcolname)