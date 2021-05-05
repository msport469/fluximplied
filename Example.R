#source('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/fluximplied.R')
inputformat='df'
species='Mmu'
geneformat='Symbol'
padjcolname='weighted_pvalue'
exampledeseqresultdataframe<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/exampledeseqresultdataframe.csv',row.names = c(1))
inputdat=exampledeseqresultdataframe

fluximplied(input,
            species,
            geneformat,
            inputformat,
            padjcolname)







#uncomment the line below this and load this for an exampledeseq result. This is the deseqresult that was put in as the argument to "as.data.frame()". 
#See code that generated this example below.
#exampledeseqresultdataframe<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/exampledeseqresultdataframe.csv',row.names = c(1))
#You would want to subset this data frame on either positive or negative LFC values, to determine whether or not the effects are positive or negative. There is a filter built into the function to subset on padjadj (adjusting the adjusted p values when we test for significance) after we perform our statistical test. 

#The code from this line up until the line below was used to generate the example
#example<-readRDS('exampledeseqresult2.RDS')
#exampledf<-as.data.frame(example)
#write.csv(exampledf,'exampledeseqresultdataframe.csv')
#------------------------------------------------------

#UNCOMMENT EVERYTHING BELOW THIS TO TEST THE FUNCTION
#input=c('Tnfa','Cpt1a')
#input=exampledf
#inputformat='vector'
#species='Mmu'
#geneformat='Symbol'
#padjcolname='weighted_pvalue'
#fluximplied(input,
#            species,
#            geneformat,
#            inputformat,
#            padjcolname)
#input=c('Tnfa','Cpt1a')
#input=exampledeseqresultdataframe
#inputup<-subset(exampledeseqresultdataframe,exampledeseqresultdataframe$log2FoldChange > 0)
#inputdown<-subset(exampledeseqresultdataframe,exampledeseqresultdataframe$log2FoldChange < 0)