#create a function that makes sure the user formatted their data according to species and format correctly. 
#If you tell the function you are giving us symbols when you are not, or otherwise lie to the function, it will give spurious results.
specform <- function(species,format) {
  ifelse(species=='Mmu', 
         ifelse(format=='Symbol', RLSgenes<-RLSdatabase$mouse.gene.symbol, 
                ifelse(format=='ENTREZID', RLSgenes<-RLSdatabase$mouse.entrez, print('Your species was accepted as Mmu, but your format was neither Symbol nor ENTREZID'))), 
         ifelse(species=='Hsa',
                ifelse(format=='Symbol', RLSgenes<-RLSdatabase$human.gene.symbol,
                       ifelse(format=='ENTREZID', RLSgenes<-RLSdatabase$human.entrez,print('Your species was accepted as Hsa, but your format was neither Symbol nor ENTREZID'))),
                print('Your species must be Mmu or Hsa')))
  RLSgenes<<-RLSgenes
}

fluximplied <- function(genes,species,format) {
  # function to see if there are any rate limiting steps in gene list
  #load the rate limiting step database
  RLSdatabase<-read.csv('https://raw.githubusercontent.com/msport469/fluximplied/master/RLSdatabase.csv',stringsAsFactors = F)
  #convert the database that matches your data for species and format (Symbol or ENTREZID)
  RLSgenes<-specform(species,format)
  #save the pathways
  RLSpathways<-RLSdatabase$metabolic.reaction
  #make a dataframe with those genes and pathways
  RLS<-data.frame(RLSgenes,RLSpathways)
  #subset the database to onlyinclude genes in your set
  subset<-subset(RLS,RLS$RLSgenes %in% genes)
  #change the column names so the user knows what each column actually is
  colnames(subset)<-c('RLS genes in your set','Pathway associated with gene')
  #create an intersect so we can actually count them
  intersect<-intersect(genes,RLSgenes)
  lengthintersect<-length(intersect)
  #create the responses if there are any genes left after subsetting on your genes
  #that are also in our database for being rate limiting steps
  ifzero<-paste('There are no genes in your set that are in our rate limiting step database. Make sure you gave the correct format of species (Mmu or Hsa only) and format (Symbol or ENTREZID only). Your genes should be in a character vector')
  ifnotzero<-paste0('Your gene set has --------> ',lengthintersect,' <-------- genes that have been identified as encoding enzymes involved as rate-limiting steps in the gene set you provided.')
  print(ifelse(lengthintersect==0,ifzero,ifnotzero))
  print(ifelse(lengthintersect==0,'Sorry about that.','Your RLS genes are saved as myRLSgenes and a dataframe of genes and corresponding pathways is saved as myRLStable'))
  #save the outputs so the user can hold onto them and look at them
  myRLStable<<-subset
  myRLSgenes<<-intersect(RLS$RLSgenes,genes)
  #print the RLS database that has been subset to only include genes that are in user's list
  subset
}

#this is an example function
genes<-c('Idh2','Cpt1a','Ifng')
species<-'Mmu'
format<-'Symbol'
fluximplied(genes=genes,species=species,format=format)
