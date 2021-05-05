#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(Cairo)
if (interactive()) {
  source("./fluximplied.R") 
# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("shinyflux"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            #fileInput("datafile", "Choose CSV File", accept = ".csv"),
            #checkboxInput("header", "Header", TRUE),
            selectInput("species", "Species", 
                                 c('Mouse'='Mmu',
                                   'Human'='Hsa')),
            selectInput('geneformat','Gene format',
                                    c('Symbol','ENTREZID')),
            selectInput('inputformat','Input format',
                                     c('Dataframe','Vector')),
            selectInput("padjcolname", "Column with adjusted p values",c(colnames(inputdat))),
           # sliderInput("Padjadj cutoff", "Padjadj cutoff",
          #                       min = 0, max = 0.2,
           #                      value = .05)
        ),

        # Show a plot of the generated distribution
        mainPanel(
          textOutput(outputId = "text"),
          textOutput(outputId = "text2"),
          
          
          tableOutput(outputId = 'table'),
          tableOutput(outputId = 'table2'),
          #tableOutput(outputId = 'significancetable')
          #textOutput(outputId = 'fluxshiny'),
          tableOutput(outputId = 'significancetable'),
          plotOutput(outputId = 'plot')
        )
    )
)

# Define server logic required to draw a histogram
#exampledeseqresultdataframe<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/exampledeseqresultdataframe.csv',row.names = c(1))
#inputdat=exampledeseqresultdataframe
server <- function(input, output) {
  output$table <- renderTable({
    data.frame(input$species,input$geneformat,input$inputformat,input$padjcolname)
  })
  #output$table2 <- renderTable(head({inputdat}),rownames = T)
  output$significancetable <- renderTable(significancetable,rownames = T)
  output$plot<-renderPlot(fluximpliedplot)
  #shinyfluximplied<-reactive({fluximplied(inputdat=inputdat,species=input$species,
  #                                        geneformat=input$geneformat,
  #                                        inputformat=input$inputformat,
  #                                        padjcolname+input$padjcolname,
  #                                        pcutoff=input$pcutoff
  #)}) 
  #output$out<-shinyfluximplied()
  

  
  #  output$text <- renderText({
  #  as.character(paste(input$pcutoff))
  #})
#  dataframe<-reactive({
#    if (is.null(input$datafile))
#      return(NULL)                
#    inputdat<-read.csv(input$datafile$datapath,row.names = c(1))
#    })
#  output$table <- renderTable({
#    head(inputdat)
#  })
  
#  

# output$answer <-renderTable({
#   shinfluximp<-shinyfluximplied()
#   #shinfluximp<-unlist(shinfluximp)
#   print(shinfluximp[[1]])})
  
 # fileInputin <-reactive({input$file1
#  })
#  
#  speciesin <- reactive({(input$species)
#  })
#  geneformatin <- reactive({(input$geneformat)
#  })
#  inputformatin <- reactive({
#    switch(input$inputformat,
#           "Dataframe" = Dataframe,
#           "Vector" = Vector)
#  })
#  padjcolin <- reactive({input$padjcol
#  })
#  padjcutoffin <- reactive({input$padjcutoff
#  })
#  
#  df <- reactive({
#    table(speciesin,geneformatin)
#  })
#  
#  output$table<-renderTable(df())
  
 # speciesinput <- reactive({
 #   switch(input$species,
 #          'Mouse'=Mmu,
 #          'Human'=Hsa)})
#  speciesinput <- reactive({
#    switch(input$species)})
  

 #   output$fluximpliedout<-fluximplied(inputdat=input$file1,
 #               species=speciesinput,
  #              geneformat=input$geneformat,
  #              inputformat=input$inputformat,
  #              padjcolname=input$padjcolname,
  #              pcutoff=input$pcutoff)
#    output$distPlot <- renderPlot({
#        ggplot(significancetable, aes(x=reorder(metabolicrxn,log2FoldChange), y=log2FoldChange , label=log2FoldChange)) + 
#            geom_bar(stat='identity', aes(fill=padjadj), width=.5,position="dodge")  +
#            scale_fill_viridis() + 
#            labs(title= "Pathway analysis with 'fluximplied'",x='Metabolic pathway',y='Log fold change',fill='Padjadj') +
#            theme(axis.title = element_text(size=12),
#                  axis.text = element_text(size=12))+
#            coord_flip()
#    })

}}
# Run the application 

shinyApp(ui = ui, server = server)


