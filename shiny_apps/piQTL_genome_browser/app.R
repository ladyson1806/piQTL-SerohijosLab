library(shiny)
library(R.utils)
library(glue)
library(igvShiny)
library(BiocManager)

options(repos = BiocManager::repositories())

PPI_table <- read.csv("./data/PPI_reference_barcodes.csv")

ui <- fluidPage(
  sidebarPanel(
    height = 6,
    column(
      6,
      selectInput(
        inputId = "ppi",
        h4("PPI of interest"),
        choices = PPI_table$PPI
      )
    ),
    column(
      6,
      selectInput(
        inputId = "drug",
        h4("Drug of interest"),
        choices = c("Fluconazole", "5-FC", "Metformin", "Trifluoperazine")
      )
    ),
    width = 12,
    h4("Loading annotations"),
    actionButton("addTracks_0", "Add LD blocks, CUT, SUT and XUT annotations"),
    br(), br(),
    actionButton("addGwasTrackButton_DrugMTX", "Add piQTL Track (DRUG)"),
    actionButton("addGwasTrackButton_noDrugMTX", "Add piQTL Track (noDRUG)"),
    br(), br(),
    h4("Region of interest"),
    textInput("roi", label = "", placeholder = "Gene or chrN:start-end"),
    actionButton("searchButton", "Search")
  ),
  
  mainPanel(
    width = 12,
    igvShinyOutput("igvShiny_0")
  )
)

server <- function(input, output, session) {
  
  annotation_base_url <- "https://raw.githubusercontent.com/ladyson1806/piQTL-SerohijosLab/main/shiny_apps/piQTL_genome_browser/data/genome_annotations"
  
  observeEvent(input$searchButton, {
    searchString <- isolate(input$roi)
    printf("--- search: %s", searchString)
    if (nchar(searchString) > 0) {
      showGenomicRegion(session, id = "igvShiny_0", searchString)
    }
  })
  
  output$igvShiny_0 <- renderIgvShiny({
    cat("--- starting renderIgvShiny\n")
    
    genomeOptions <- parseAndValidateGenomeSpec(genomeName = "sacCer3")
    
    x <- igvShiny(
      genomeOptions,
      tracks = list()
    )
    
    cat("--- ending renderIgvShiny\n")
    x
  })
  
  observeEvent(input$addTracks_0, {
    printf("---- add remote GFF3 tracks")
    
    loadGFF3TrackFromURL(
      session, id = "igvShiny_0",
      trackName = "LD blocks",
      gff3URL = glue("{annotation_base_url}/LD_blocks_050.gff3.gz"),
      indexURL = glue("{annotation_base_url}/LD_blocks_050.gff3.gz.tbi"),
      color = "grey",
      colorByAttribute = "type",
      colorTable = list(),
      displayMode = "EXPANDED",
      trackHeight = 50,
      visibilityWindow = 100000
    )
    
    loadGFF3TrackFromURL(
      session, id = "igvShiny_0",
      trackName = "CUTs Xu 2009",
      gff3URL = glue("{annotation_base_url}/Xu_2009_CUTs_V64.gff3.gz"),
      indexURL = glue("{annotation_base_url}/Xu_2009_CUTs_V64.gff3.gz.tbi"),
      color = "darkgreen",
      colorByAttribute = "type",
      colorTable = list(),
      displayMode = "EXPANDED",
      trackHeight = 50,
      visibilityWindow = 100000
    )
    
    loadGFF3TrackFromURL(
      session, id = "igvShiny_0",
      trackName = "SUTs Xu 2009",
      gff3URL = glue("{annotation_base_url}/Xu_2009_SUTs_V64.gff3.gz"),
      indexURL = glue("{annotation_base_url}/Xu_2009_SUTs_V64.gff3.gz.tbi"),
      color = "green",
      colorByAttribute = "type",
      colorTable = list(),
      displayMode = "EXPANDED",
      trackHeight = 50,
      visibilityWindow = 100000
    )
    
    loadGFF3TrackFromURL(
      session, id = "igvShiny_0",
      trackName = "XUTs Van Dijk 2011",
      gff3URL = glue("{annotation_base_url}/van_Dijk_2011_XUTs_V64.gff3.gz"),
      indexURL = glue("{annotation_base_url}/van_Dijk_2011_XUTs_V64.gff3.gz.tbi"),
      color = "palegreen",
      colorByAttribute = "type",
      colorTable = list(),
      displayMode = "EXPANDED",
      trackHeight = 50,
      visibilityWindow = 100000
    )
  })

  observeEvent(input$addGwasTrackButton_DrugMTX, {
    track_url <- as.character(glue(
      "https://raw.githubusercontent.com/ladyson1806/public_hosting/main/piQTL_mapping/formatted_for_genome_browser/{input$ppi}_{input$drug}_SNP_annotations.tsv"
    ))
    
    track <- GWASTrack(
      glue("{input$ppi} with {input$drug}"),
      track_url,
      chrom.col = 2,
      pos.col = 3,
      pval.col = 5,
      trackHeight = 100
    )
    
    display(track, session, id = "igvShiny_0", deleteTracksOfSameName = TRUE)
  })
  
  observeEvent(input$addGwasTrackButton_noDrugMTX, {
    track_url <- as.character(glue(
      "https://raw.githubusercontent.com/ladyson1806/public_hosting/main/piQTL_mapping/formatted_for_genome_browser/{input$ppi}_noDrug_SNP_annotations.tsv"
    ))
    
    track <- GWASTrack(
      glue("{input$ppi} without {input$drug}"),
      track_url,
      chrom.col = 2,
      pos.col = 3,
      pval.col = 5,
      trackHeight = 100
    )
    
    display(track, session, id = "igvShiny_0", deleteTracksOfSameName = TRUE)
  })
  }

shinyApp(ui = ui, server = server)