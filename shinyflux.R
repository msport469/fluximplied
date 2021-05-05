#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            inputdat<-inputdat,
            species<-selectInput("species", "Species", 
                                 c('Mouse'='Mmu',
                                   'Human'='Hsa')),
            geneformat<-selectInput('geneformat','Gene format',
                                    c('Symbol','ENTREZID')),
            inputformat<-selectInput('inputformat','Input format',
                                     c('Vector','Dataframe')),
            padjcolname<-selectInput("padjcolname", "Column with adjusted p values",c(colnames(inputdat))),
            pcutoff<-sliderInput("Padjadj cutoff", "Padjadj cutoff",
                                 min = 0, max = 0.2,
                                 value = .05)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    fluximplied(inputdat=input$inputdat,
                species,
                geneformat,
                inputformat,
                padjcolname,
                pcutoff)
#    output$distPlot <- renderPlot({
#        ggplot(significancetable, aes(x=reorder(metabolicrxn,log2FoldChange), y=log2FoldChange , label=log2FoldChange)) + 
#            geom_bar(stat='identity', aes(fill=padjadj), width=.5,position="dodge")  +
#            scale_fill_viridis() + 
#            labs(title= "Pathway analysis with 'fluximplied'",x='Metabolic pathway',y='Log fold change',fill='Padjadj') +
#            theme(axis.title = element_text(size=12),
#                  axis.text = element_text(size=12))+
#            coord_flip()
#    })
}

# Run the application 
shinyApp(ui = ui, server = server)
