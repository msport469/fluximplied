dp49a<-read.csv('deseq2result.dpvscd49a.csv') #,row.names = c(1)
head(dp49a)
dim(dp49a)
updp49a<-subset(dp49a,dp49a$log2FoldChange>0)
head(updp49a)
dim(updp49a)
downdp49a<-subset(dp49a,dp49a$log2FoldChange<0)
head(downdp49a)
dim(downdp49a)



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









#input=downdp49a
inputformat='df'
species='Mmu'
geneformat='Symbol'
padjcolname='weighted_pvalue'









#create a function that makes sure the user formatted their data according to species and format correctly. 
#If you tell the function you are giving us symbols when you are not, or otherwise lie to the function, it will give spurious results.
specform <- function(species,geneformat) {
  RLSdatabase<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/RLSdatabase.csv',stringsAsFactors = F)
  ifelse(species=='Mmu' || species=='MMU' || species=='mmu', 
         ifelse(geneformat=='Symbol', RLSgenes<-RLSdatabase$mouse.gene.symbol, 
                ifelse(geneformat== 'Entrezid'||geneformat=='ENTREZ'||geneformat=='Entrez'||geneformat=='entrez'||geneformat=='entrezid'||geneformat=='ENTREZID',
                       RLSgenes<-RLSdatabase$mouse.entrez, print('Your species was accepted as Mmu, but your geneformat was neither Symbol nor ENTREZID'))), 
         ifelse(species=='Hsa'||species=='hsa'||species=='HSA',
                ifelse(geneformat=='Symbol'||geneformat=='SYMBOL'||geneformat=='symbol', RLSgenes<-RLSdatabase$human.gene.symbol,
                       ifelse(geneformat== 'Entrezid'||geneformat=='ENTREZ'||geneformat=='Entrez'||geneformat=='entrez'||geneformat=='entrezid'||geneformat=='ENTREZID',
                              RLSgenes<-RLSdatabase$human.entrez,print('Your species was accepted as Hsa, but your geneformat was neither Symbol nor ENTREZID'))),
                print('Your species must be Mmu or Hsa')))
  RLSgenes<<-RLSgenes
}

fluximplied <- function(input,species,geneformat,inputformat,padjcolname) {
  # function to see if there are any rate limiting steps in gene list
  #load the rate limiting step database
  RLSdatabase<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/RLSdatabase.csv',stringsAsFactors = F)
  #convert the database that matches your data for species and geneformat (Symbol or ENTREZID)
  RLSgenes<-specform(species,geneformat)
  #save the pathways
  RLSpathways<-RLSdatabase$metabolic.reaction
  #make a dataframe with those genes and pathways
  RLS<-data.frame(RLSgenes,RLSpathways)
  #create the responses if there are any genes left after subsetting on your genes
  #that are also in our database for being rate limiting steps
   ifelse(inputformat=='vector'||inputformat=='Vector'||inputformat=='VECTOR',
         print('We are using your vector of genes as the input'),
         ifelse(inputformat=='df'||inputformat=='DF'||inputformat=='Df'||inputformat=='dataframe'||inputformat=='Dataframe',{
                input<-subset(input,rownames(input) %in% RLSgenes) 
                  input$padjadj<-p.adjust(input[[padjcolname]],method = 'BH') 
                  inputssubset<-subset(input,input$padjadj<.05) 
                  input<-rownames(inputssubset)},
                print('It appears that you supplied an input that was neither a dataframe nor a vector.')
         )
  )
  #subset the database to only include genes in your set
  subset<-subset(RLS,RLS$RLSgenes %in% input)
  #change the column names so the user knows what each column actually is
  colnames(subset)<-c('RLS genes in your set','Pathway associated with gene')
  #create an intersect so we can actually count them
  intersect<-intersect(input,RLSgenes)
  lengthintersect<-length(intersect)
  print(ifelse(lengthintersect==0,
               paste('There are no genes in your set that are in our rate limiting step database. Make sure you gave the correct species (Mmu or Hsa only) and geneformat (Symbol or ENTREZID only). Your genes should be in a character vector. Sorry about that. We are as sad as you.'),
               {paste0('Your gene set has --------> ',lengthintersect,' <-------- genes that have been identified as encoding enzymes involved as rate-limiting steps in the gene set you provided. Your RLS genes are saved as myRLSgenes and a dataframe of genes and corresponding pathways is saved as myRLStable')}))
  #save the outputs so the user can hold onto them and look at them
<<<<<<< HEAD
  myRLStable<<-subset
  myRLSgenes<<-intersect(RLS$RLSgenes,input)
  #print the RLS database that has been subset to only include genes that are in user's list
  subset
=======
   myRLStable<<-subset
   myRLSgenes<<-intersect(RLS$RLSgenes,input)
  ifelse(inputformat=='df'||inputformat=='DF'||inputformat=='Df'||inputformat=='dataframe'||inputformat=='Dataframe',
         {significancetable<-inputssubset
         significancetable$metabolicrxn <- myRLStable$`Pathway associated with gene`[match(rownames(significancetable), myRLStable$`RLS genes in your set`)]
         significancetable<<-significancetable
         ifelse (!require(ggplot2),stop("ggplot2 not installed"),1+1)
         ifelse (!require(viridis),stop("viridis not installed"),1+1)
         fluximpliedplot<<-ggplot(significancetable, aes(x=reorder(metabolicrxn,log2FoldChange), y=log2FoldChange , label=log2FoldChange)) + 
           geom_bar(stat='identity', aes(fill=padjadj), width=.5,position="dodge")  +
           scale_fill_viridis() + 
           labs(title= "Pathway analysis with 'fluximplied'",x='Metabolic pathway',y='Log fold change',fill='Padjadj') +
           theme(axis.title = element_text(size=12),
                 axis.text = element_text(size=12))+
           coord_flip()
         plot(fluximpliedplot)},
         print(subset))

>>>>>>> 2fd4da13ef2b6d1a778ed740ca4b167a114f61d7
}


input=updp49a
inputformat='df'
species='Mmu'
geneformat='Symbol'
padjcolname='weighted_pvalue'

fluximplied(input,
            species,
            geneformat,
            inputformat,
            padjcolname)



names(dp49a)[1] <- "genes"
head(dp49a)
sub<-subset(dp49a, dp49a$genes %in% RLSdatabase$mouse.gene.symbol)
sub


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
<<<<<<< HEAD
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
=======
input=c('Tnfa','Cpt1a')
input=exampledeseqresultdataframe
inputup<-subset(exampledeseqresultdataframe,exampledeseqresultdataframe$log2FoldChange > 0)
inputdown<-subset(exampledeseqresultdataframe,exampledeseqresultdataframe$log2FoldChange < 0)

inputformat='df'
species='Mmu'
geneformat='Symbol'
padjcolname='weighted_pvalue'

fluximplied(input,
            species,
            geneformat,
            inputformat,
            padjcolname)

>>>>>>> 2fd4da13ef2b6d1a778ed740ca4b167a114f61d7
