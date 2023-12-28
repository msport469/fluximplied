# Install and load package
install.packages('devtools')
library(devtools)
install_github('sportiellomike/fluximplied',build_vignettes=T)
library(fluximplied)

# data frame like input
inputdat<-exampleData
fluximplied(inputdat = inputdat,
            species = 'mmu',
            geneformat = 'symbol',
            inputformat = 'dataframe',
            padjcolname = 'padj',
            pcutoff = .05
              )
ggsave('./exampleplot.png',plot = fluximpliedplot,device = 'tiff',dpi = 600,units = 'in',height = 4,width=8)
#vector like input
inputdat<-c('Ifng','Idh2','Pfkl')
fluximplied(inputdat = inputdat,
            species = 'mmu',
            geneformat = 'symbol',
            inputformat = 'vector'
            )
