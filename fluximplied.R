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

fluximplied <- function(inputdat,species,geneformat,inputformat,padjcolname='adj_pvalue',pcutoff=0.05) {
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
         print('We are using your vector of genes as the inputdat'),
         ifelse(inputformat=='df'||inputformat=='DF'||inputformat=='Df'||inputformat=='dataframe'||inputformat=='Dataframe',{
                inputdat<-subset(inputdat,rownames(inputdat) %in% RLSgenes) 
                  inputdat$padjadj<-p.adjust(inputdat[[padjcolname]],method = 'BH') 
                  inputssubset<-subset(inputdat,inputdat$padjadj<pcutoff) 
                  inputdat<-rownames(inputssubset)},
                print('It appears that you supplied an input that was neither a dataframe nor a vector.')
         )
  )
  #subset the database to only include genes in your set
  subset<-subset(RLS,RLS$RLSgenes %in% inputdat)
  #change the column names so the user knows what each column actually is
  colnames(subset)<-c('RLS genes in your set','Pathway associated with gene')
  #create an intersect so we can actually count them
  intersect<-intersect(inputdat,RLSgenes)
  lengthintersect<-length(intersect)
  print(ifelse(lengthintersect==0,
               paste('There are no genes in your set that are in our rate limiting step database. Make sure you gave the correct species (Mmu or Hsa only) and geneformat (Symbol or ENTREZID only). Your genes should be in a character vector. Sorry about that. We are as sad as you.'),
               {paste0('Your gene set has --------> ',lengthintersect,' <-------- genes that have been identified as encoding enzymes involved as rate-limiting steps in the gene set you provided. Your RLS genes are saved as myRLSgenes and a dataframe of genes and corresponding pathways is saved as myRLStable')}))
  #save the outputs so the user can hold onto them and look at them
  myRLStable<<-subset
  myRLSgenes<<-intersect(RLS$RLSgenes,inputdat)
  #print the RLS database that has been subset to only include genes that are in user's list
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
         plot(fluximpliedplot)}#,
        # print(subset)
        )
}

