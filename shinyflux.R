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
fields <- c('inputdat',
            "species",
            "geneformat",
            "inputformat",
            'padjcolname',
            'pcutoff')
#exampledeseqresultdataframe<-read.csv('https://raw.githubusercontent.com/sportiellomike/fluximplied/master/exampledeseqresultdataframe.csv',row.names = c(1))
inputdat=exampledeseqresultdataframe
source("./fluximplied.R")
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
ui <- fluidPage(
  shinythemes::themeSelector(),  # <--- Add this somewhere in the UI

    # Application title
    titlePanel("shinyfluximplied"),
    
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            #fileInput("datafile", "Choose CSV File", accept = ".csv"),
            #checkboxInput("header", "Header", TRUE),
          selectInput("inputdat", "inputdat", 
                      c('inputdat'='inputdat')),  
          selectInput("species", "species", 
                                 c('Mouse'='Mmu',
                                   'Human'='Hsa')),
            selectInput('geneformat','Gene format',
                                    c('Symbol','ENTREZID')),
            selectInput('inputformat','Input format',
                                     c('Dataframe','Vector')),
            selectInput("padjcolname", "Column with adjusted p values",c(colnames(inputdat))),
            sliderInput("pcutoff", "Padjadj cutoff",
                                 min = 0, max = 0.2,
                                 value = .05),
            actionButton("submit", "Submit")
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
  
output$table <- renderTable({
  species <-input$species
  geneformat <-input$geneformat
  inputformat <-input$inputformat
  padjcolname <-input$padjcolname
  pcutoff <- input$pcutoff
  fluximplied(inputdat,species,geneformat,inputformat,padjcolname,pcutoff)
  significancetable
}rownames = T)
output$plot <-renderPlot({
  species <-input$species
  geneformat <-input$geneformat
  inputformat <-input$inputformat
  padjcolname <-input$padjcolname
  pcutoff <- input$pcutoff
  fluximplied(inputdat,species,geneformat,inputformat,padjcolname,pcutoff)
  fluximpliedplot
})
output$print <-renderText({
  species <-input$species
  geneformat <-input$geneformat
  inputformat <-input$inputformat
  padjcolname <-input$padjcolname
  pcutoff <- input$pcutoff
  fluximplied(inputdat,species,geneformat,inputformat,padjcolname,pcutoff)
})
 
}}

shinyApp(ui = ui, server = server)


