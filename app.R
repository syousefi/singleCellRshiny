library(shiny)
library(DT)
library(Seurat)
source("global.R")

options(shiny.maxRequestSize= 180*1024^2) # defult upload size is 5MB, I increased to 180 MB

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Single-cell RNA-seq Data Analysis!"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Sidebar layout for inputing csv files ---
      # Sidebar panel for inputs ----
      fileInput('file1', h5('Upload your single-cell dataset', alighn= "left"), accept= c(".csv") ), 
      textInput("gene", h5("Select the gene of interest", align = "left"), value= "Thy1") ,
      
      sliderInput(inputId = "minCells",
                  label = "Min number of cells expressed a particular gene:",
                  min = 1,
                  max = 50,
                  value = 30),
      
      sliderInput(inputId = "minGenes",
                  label = "Min number of genes expressed in one cell:",
                  min = 100,
                  max = 1000,
                  value = 200),
      
      br(),
      
      sliderInput(inputId = "xLowCutOff",
                  label = "Average expression low cut-off:",
                  min = 0,
                  max = 2, step = 0.1, 
                  value = 0.1),
      
      br(),
      
      sliderInput(inputId = "yLowCutOff",
                  label = "Dispersion low cut-off:",
                  min = -4,
                  max = 5,
                  value = 2),
      
      
      br(),
      
      actionButton("nGenesHist", "Histogram of the number of genes"),
      
      br(),
      
      actionButton("rgcHist", "Histogram of the number of RGCs"),
      
      br(),
      
      actionButton("expDispPlot", "Average expression versus dispersion"),
      
      br(),
      
      # actionButton("update", "Update view")
      h5("OPAR single-cell data analysis"),
      
      br(),  # put a line space
      
      submitButton("Run")
      
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram of number of genes ----
      # plotOutput(outputId = "genePlot")
      
      # DT::dataTableOutput('contents'),
      # br(),
      tabsetPanel(
        tabPanel("Insight",
                 # Text: filtering genes and cells
                 h5(textOutput("insight"),
                    
                    tableOutput("view"),
                    
                    verbatimTextOutput("summary")),
                 
                 # Text: average gene expression text output
                 h5(textOutput("aveGene")) ),
        
        tabPanel("Number of Genes",
                 # Histogram: number of genes histogram
                 plotOutput('nGenePlot')),
        
        tabPanel("RGC markers", 
                 # Histogram: number of RGC markers histogram
                 plotOutput("nRGCmarkersPlot")),
        tabPanel("Variable genes",
                 # Scatter plot: average expression versus dispersion plot
                 plotOutput("expDispPlot", width= 1600, height = 900))
        
      )
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # data <- reactiveValues()
  # # to read the csv dataset
  # inputData <- renderTable({
  #   
  #   # read the dataset
  #   req(input$file1)
  #   inFile <- input$file1
  #   read.csv(inFile$datapath, head= TRUE, sep= ",", stringsAsFactors= FALSE)
  #   assign('data', data2, envir=.GlobalEnv)
  #   print(summary(data))
  #   
  # })
  
  inputData <- reactive({
    inFile <- input$file1
    req(input$file1)
    data <- read.csv(inFile$datapath, head= TRUE, sep= ",", stringsAsFactors= FALSE)
    return(data)
  })
  
  output$contents <- DT::renderDataTable({
    DT::datatable(inputData())
    
  })
  
  # print some dataset info
  output$summary <- renderPrint({
#    nrow(inputData())})
   summary(inputData()$H_1_S41_ROW01)})
# 
# 
  # Show the first "n" observations
  output$view <- renderTable({
    head(inputData())
  })

  # this will print a text
  output$insight <- renderText({
    # "You have selected this"
    sc.object.fly  <- CreateSeuratObject(raw.data= data.log, min.cells= input$minCells, min.genes= input$minGenes, project= "Single_cell_data" )  # 2000 genes, 900 genes in Macosko, genes expressed in fewer than 0.05% of cells are exluded
    paste0(toString(dim(data.log)[2]-dim(sc.object.fly@raw.data)[2]) , " cells (out of 800) and ", toString(dim(data.log)[1]-dim(sc.object.fly@raw.data)[1]), " genes (out of 25,394) were excluded")
  })
  
  output$aveGene <- renderText({
    aveGeneExp= paste0("Average expression of the selected gene is: ", ave.expressed.genes.2[input$gene])
    
  })
#   
  
  # to plot the histogram of the number of genes per cell 
  observeEvent(input$nGenesHist, {   # this is to make this plot dependent on this button
    output$nGenePlot <- renderPlot({
      hist(n.genes.per.cell, breaks= 100, col= c("gray"), 
           xlab = "Number of genes expressed", # main= "Histogram of number of genes expressed per cell",
           cex.axis= 1.5, cex.lab= 1.6)
    })
  })
  
  # to plot the histogram of the number of RGC maekers per cell
  observeEvent(input$rgcHist, {   # this is to make this plot dependent on this button
    output$nRGCmarkersPlot <- renderPlot({
      barplot(can.cell.counts, main="", xlab="Number of known RGC genes expressed", ylab= "", col= c("gray"), 
              cex.axis= 1.5, cex.names=1.5, cex.lab= 1.8 ) # title: Cells with different number of reference RGC genes
      title(ylab= "Frequencey", line= 2.5, cex.lab= 1.8)
      grid()
    })
  })
  
  # to plot the expression of genes versus disperssion
  observeEvent(input$expDispPlot, {
    output$expDispPlot <- renderPlot({
      sc.object <- FindVariableGenes(object= sc.object, mean.function= ExpMean, dispersion.function= LogVMR, x.low.cutoff= input$xLowCutOff, x.high.cutoff= Inf, y.cutoff= input$yLowCutOff, num.bin= 20, binning.method= "equal_width", do.plot= FALSE)
      # plot( sc.object@hvg.info$gene.mean, sc.object@hvg.info$gene.dispersion.scaled )
      # average
      # plot(input$xLowCutOff, input$yLowCutOff)
      VariableGenePlot(sc.object, do.text= TRUE, cex.use= 0.5, cex.text.use= 0.9, do.spike= FALSE, pch.use= 19, col.use= "gray",  spike.col.use= "red", plot.both= FALSE, do.contour= TRUE,
                        contour.lwd= 3, contour.col= "white", contour.lty= 2, x.low.cutoff= input$xLowCutOff, x.high.cutoff= Inf, y.cutoff= input$yLowCutOff)
      pass.cutoff.1= which( rownames(sc.object@hvg.info) %in% input$gene ) # which( rownames(sc.object@hvg.info) %in% c(reference.rgc.markers.4, subtypes.rgc.markers.4) ) # Rheaume pan and subtype RGCs
      pass.cutoff.2= which( rownames(sc.object@hvg.info) %in% c(reference.rgc.markers.3) )                         # Rob's RGC genes
      text( sc.object@hvg.info$gene.mean[pass.cutoff.1], sc.object@hvg.info$gene.dispersion.scaled[pass.cutoff.1], rownames(sc.object@hvg.info)[pass.cutoff.1], col= "red", cex= 2.0)
      text( sc.object@hvg.info$gene.mean[pass.cutoff.2], sc.object@hvg.info$gene.dispersion.scaled[pass.cutoff.2], rownames(sc.object@hvg.info)[pass.cutoff.2], col= "blue", cex= 0.8)
      
      
    })
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
