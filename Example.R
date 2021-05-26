#load fluximplied from github
source('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/fluximplied.R')
#load example inputdata in the format of dataframe
exampledeseqresultdataframe<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/exampledeseqresultdataframe.csv',row.names = c(1))
inputdat=exampledeseqresultdataframe
#below is an example input dat that is a character vector
#inputdat=c('Tnfa','Cpt1a')
#load other parameters
inputformat='df'
species="mmu"
geneformat="SYMBOL"
padjcolname='weighted_pvalue'
pcutoff=0.05
#now actually run fluximplied
fluximplied(inputdat,
            species,
            geneformat,
            inputformat,
            padjcolname,
            pcutoff)

#download pathways with pathview package.

###FIN###

##Other stuff##
#uncomment the line below this and load this for an exampledeseq result. This is the deseqresult that was put in as the argument to "as.data.frame()". 
#See code that generated this example below.
#exampledeseqresultdataframe<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/exampledeseqresultdataframe.csv',row.names = c(1))
#You would want to subset this data frame on either positive or negative LFC values, to determine whether or not the effects are positive or negative. There is a filter built into the function to subset on padjadj (adjusting the adjusted p values when we test for significance) after we perform our statistical test. 

#The code from this line was used to generate the example CSV
#example<-readRDS('exampledeseqresult2.RDS')
#exampledf<-as.data.frame(example)
#write.csv(exampledf,'exampledeseqresultdataframe.csv')