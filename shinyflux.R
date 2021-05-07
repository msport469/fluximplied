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
library(shinythemes)
#fields <- c('inputdat',
#            "species",
#            "geneformat",
#            "inputformat",
#            'padjcolname',
#            'pcutoff')
#exampledeseqresultdataframe<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/exampledeseqresultdataframe.csv',row.names = c(1))
#inputdat=exampledeseqresultdataframe
source("https://raw.githubusercontent.com/sportiellomike/fluximplied/master/fluximplied.R")
saveData <- function(data) {
  data <- as.data.frame(t(data))
  if (exists("responses")) {
    responses <<- rbind(responses, data)
  } else {
    responses <<- data
  }
}

loadData <- function() {
  if (exists("responses")) {
    responses
  }
}
if (interactive()) {
# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("slate"),
    # Application title
    titlePanel("fluximplied"),
    
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("file1", "Choose CSV File", accept = ".csv"),
            #checkboxInput("header", "Header", TRUE),
         # selectInput("inputdat", "inputdat", 
          #            c('inputdat'='inputdat')),  
          selectInput("species", "species", 
                                 c('Mouse'='Mmu',
                                   'Human'='Hsa')),
            selectInput('geneformat','Gene format',
                                    c('Symbol','ENTREZID')),
            selectInput('inputformat','Input format',
                                     c('Dataframe','Vector')),
            selectInput("padjcolname", "Column with adjusted p values",''),
            numericInput("pcutoff", "Significance cutoff (alpha)", 0.05, min = 0, max = 1),
            downloadButton("downloadData", "Download output table")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tableOutput(outputId = 'table'),
          plotOutput(outputId = 'plot'),
          textOutput(outputId = "print"),
        )
    )
)
server = function(input, output, session) {
#  outVar = reactive({
#    file <- input$file1
#    ext <- tools::file_ext(file$datapath)
#    req(file)
#    validate(need(ext == "csv", "Please upload a csv file"))
#    inputdat<-read.csv(file$datapath, row.names=c(1))
#    mydata = get(inputdat)
#    names(mydata)
#  })
#  observe({
#    updateSelectInput(session, "padjcolname",
#                      choices = outVar()
#    )})
  data <- reactive({ 
    req(input$file1) ## ?req #  require that the input is available
    
    inFile <- input$file1 
    df <- read.csv(inFile$datapath, header = T,row.names = c(1))
    updateSelectInput(session, inputId = 'padjcolname', label = 'padjcol',choices = colnames(df))
    return(df)
  })
  #inputdat <- reactive({
  #  req(input$file)
  #  
  #  ext <- tools::file_ext(input$file$name)
  #  switch(ext,
  #         csv = vroom::vroom(input$file$datapath, delim = ","),
  #         tsv = vroom::vroom(input$file$datapath, delim = "\t"),
  #         validate("Invalid file; Please upload a .csv or .tsv file")
  #  )
   # read.csv(input$file$datapath, row.names=c(1))
 # })
output$table <- renderTable({
  inputdat<-data()
  species <-input$species
  geneformat <-input$geneformat
  inputformat <-input$inputformat
  padjcolname <-input$padjcolname
  pcutoff <- input$pcutoff
  fluximplied(inputdat,species,geneformat,inputformat,padjcolname,pcutoff)
  return(significancetable)
  },rownames = T)
output$plot <-renderPlot({
  inputdat<-data()
  species <-input$species
  geneformat <-input$geneformat
  inputformat <-input$inputformat
  padjcolname <-input$padjcolname
  pcutoff <- input$pcutoff
  fluximplied(inputdat,species,geneformat,inputformat,padjcolname,pcutoff)
  return(fluximpliedplot)
})
output$print <-renderText({
  inputdat<-data()
  species <-input$species
  geneformat <-input$geneformat
  inputformat <-input$inputformat
  padjcolname <-input$padjcolname
  pcutoff <- input$pcutoff
  fluximplied(inputdat,species,geneformat,inputformat,padjcolname,pcutoff)
  return(print1)
})

tabledownload <- reactive({
  inputdat<-data()
  species <-input$species
  geneformat <-input$geneformat
  inputformat <-input$inputformat
  padjcolname <-input$padjcolname
  pcutoff <- input$pcutoff
  fluximplied(inputdat,species,geneformat,inputformat,padjcolname,pcutoff)
  return(significancetable)
})

output$tabledownload <- renderTable({
 
  tabledownload()
})

output$downloadData  <- downloadHandler(
  filename = function() {
    paste('significancetable', ".csv", sep = "")
  },
  content = function(file) {
    write.csv(tabledownload(), file, row.names = T)
  }
)
 
}}

shinyApp(ui = ui, server = server)


