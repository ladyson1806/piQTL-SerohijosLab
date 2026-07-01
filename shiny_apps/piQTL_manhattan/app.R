library(shiny)
library(plotly)
options(shiny.host = '0.0.0.0')
options(shiny.port = 3838)

# Define UI for app that draws a histogram ----

PPI_table = read.csv('./data/PPI_reference_barcodes.csv')

ui <- fluidPage(
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(width=12,
    # Input: Text box for side playing PPI and Drug conditions ----
      selectInput(inputId = "ppi", h4("PPI of interest"), 
                choices = PPI_table$PPI),
      selectInput(inputId = "drug", h4("Drug of interest"), 
                choices = c("Fluconazole", "5-FC", "Metformin", "Trifluoperazine"))
  ),
    
    # Main panel for displaying outputs ----
    mainPanel(width=12,
      h3("Drug"),
      br(),
      br(),
      column(
        width = 8,
        plotlyOutput(outputId = "Drug_MTX_manhattanplot")
        ),
      
      column(
        width = 4,
        plotlyOutput(outputId = "Drug_MTX_qqplot")
      ),
      
      br(),
      br(),

      h3("noDrug"),
      column(
        width = 8,
        plotlyOutput(outputId = "noDrug_MTX_manhattanplot")
      ),
      
      column(
        width = 4,
        plotlyOutput(outputId = "noDrug_MTX_qqplot")
      )
    )
  )
)

# Define server logic required to draw a histogram ----
library(glue)
library(manhattanly)

server <- function(input, output) {
  
  Drug_MTX_datasetInput <- reactive({
    PPI <- input$ppi
    DRUG <- input$drug
    infile = glue('https://raw.githubusercontent.com/ladyson1806/public_hosting/main/piQTL_mapping/formatted_for_manhattan/{PPI}_{DRUG}_piQTLs.csv')
    read.csv(infile)
  })

  noDrug_MTX_datasetInput <- reactive({
    PPI <- input$ppi
    DRUG <- input$drug
    infile = glue('https://raw.githubusercontent.com/ladyson1806/public_hosting/main/piQTL_mapping/formatted_for_manhattan/{PPI}_noDrug_piQTLs.csv')
    
    read.csv(infile)
  })
  
  
  output$Drug_MTX_manhattanplot <- renderPlotly({
    manhattanly(Drug_MTX_datasetInput(), chr='CHR', snp="SNP", gene="GENE", labelChr=c(1:16,"MT"), suggestiveline = -log10(1e-04), genomewideline=F, title=glue('{input$ppi} with {input$drug}' ))
  })
  
  output$Drug_MTX_qqplot <- renderPlotly({
    qqly(Drug_MTX_datasetInput(), chr='CHR', snp="SNP", gene="GENE")
  })
  
  output$noDrug_MTX_manhattanplot <- renderPlotly({
    manhattanly(noDrug_MTX_datasetInput(), chr='CHR', snp="SNP", gene="GENE", labelChr=c(1:16,"MT"), suggestiveline = -log10(1e-04), genomewideline=F, title=glue('{input$ppi} without {input$drug}'))
  })
  
  output$noDrug_MTX_qqplot <- renderPlotly({
    qqly(noDrug_MTX_datasetInput(), chr='CHR', snp="SNP", gene="GENE")
  })

}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
