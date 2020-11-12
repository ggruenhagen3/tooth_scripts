##########
# TJ APP #
##########
options(repos = BiocManager::repositories())
library("shiny")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("ggplot2")
# options(warn=-1)
# tj <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/tj_shiny/data/tj.rds")
tj <- readRDS("./data/tj.rds")
gene_names <- rownames(tj@assays$RNA)

ui <- fluidPage(
  
  # App title ----
  titlePanel("Paint Expression of Cells"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      selectizeInput(inputId = "gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 10)),
      textInput(inputId = "cluster", label = "Cluster to Split", value = "0", width = NULL),
      textInput(inputId = "resolution", label = "Resolution", value = "0.2", width = NULL)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      tabsetPanel(id = "tabs", type = "tabs",
                  tabPanel("UMAP", value="umap_plot", plotOutput("umap_plot", width = "100%", height="500px")),
                  tabPanel("BarPlot Summary", value="barplot", plotOutput("barplot", width = "100%", height="500px")),
                  tabPanel("BarPlot %", value="barplotAvg", plotOutput("barplotAvg", width = "100%", height="500px")),
                  tabPanel("Summary", value="summary", verbatimTextOutput("summary")),
                  tabPanel("Summary %", value="summaryAvg", verbatimTextOutput("summaryAvg"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  updateSelectizeInput(session, "gene", choices = gene_names, server = TRUE)
  
  # Helper Functions
  findOvlPCells <- function(genes) {
    clean_genes <- c()
    pos_cells <- c()
    for (gene in genes) {
      gene_lower <- tolower(gene)
      gene_upper <- toupper(gene)
      gene_title <- str_to_title(gene)
      if (gene_lower %in% gene_names) {
        gene <- gene_lower
      } else if (gene_upper %in% gene_names) {
        gene <- gene_upper
      } else if (gene_title %in% gene_names) {
        gene <- gene_title
      }
      clean_genes <- c(gene, clean_genes)
      expr <- FetchData(object = tj, vars = gene)
      pos_cells <- c(colnames(tj[, which(x = expr > 1)]), pos_cells)
    }
    
    # Find the overlap of the positive cells
    num_genes <- length(genes)
    counts <- table(pos_cells)
    ovlp_cells <- rownames(counts[counts >= num_genes])
    tj <- SetIdent(tj, cells=ovlp_cells, value="overlapping_cells")
    tj <- SetIdent(tj, cells=setdiff(WhichCells(tj), ovlp_cells), value="non_overlapping_cells")
    tj$ovlp <- tj@active.ident
    return(tj$ovlp)
  }
  
  createBarPlot <- function(average) {
    tj$ovlp <- findOvlPCells(input$gene)
    Idents(tj) <- tj$ovlp
    # tj$sample <- rep("TJ", ncol(tj))
    test <- tj[, WhichCells(tj, idents = "overlapping_cells")]
    num_clusters <- as.numeric(tail(levels(tj@meta.data$seurat_clusters), n=1))
    total_tj <- c()
    Idents(test) <- test$sample
    Idents(tj) <- tj$sample
    test$cond.cluster <- paste(Idents(test), test$seurat_clusters, sep = "_")
    tj$cond.cluster <- paste(Idents(tj), tj$seurat_clusters, sep = "_")
    Idents(tj) <- tj$cond.cluster
    Idents(test) <- test$cond.cluster
    for (i in 0:num_clusters) {
      tj_cells_in_cluster <- 0
      try(tj_cells_in_cluster <- length(WhichCells(test, idents = paste("TJ", i, sep="_"))), silent=TRUE)
      all_tj_cells_in_cluster <- length(WhichCells(tj, idents = paste("TJ", i, sep="_")))
      
      if (average == TRUE) {
        total_tj <- c(total_tj, round( (tj_cells_in_cluster/all_tj_cells_in_cluster) * 100, 2))
      } else {
        total_tj <- c(total_tj, tj_cells_in_cluster)
      }
      
    }
    df <- data.frame(condition <- c(rep("TJ", length(total_tj))),cluster_num <- c(0:num_clusters),
                     value <- c(total_tj))
    colnames(df) <- c("condition", "cluster_num", "value")
    my_title <- paste("Number of Cells Expressing", paste(input$gene, collapse = ' and '), "per Cluster")
    my_ylab <- "Number of Cells"
    if (average == TRUE) {
      my_title <- paste("% Cells Expressing", paste(input$gene, collapse = ' and '), "per Cluster")
      my_ylab <- "% Cells"
    }
    p <- ggplot(df, aes(fill=condition, x=cluster_num, y=value)) +
      geom_bar(position="dodge", stat="identity") +
      theme_minimal() +
      ggtitle(my_title) +
      xlab("Cluster") +
      ylab(my_ylab) +
      scale_x_continuous(breaks = 0:40) +
      geom_text(aes(label=value), vjust=1.6, color="black", position = position_dodge(0.9), size=3.5)
    theme_minimal()
    p
  }
  
  createOvlpPlot <- function() {
    if (length(input$gene) > 1) {
      tj$ovlp <- findOvlPCells(input$gene)
      Idents(tj) <- tj$ovlp
      DimPlot(tj, reduction="umap", group.by = "ident", pt.size=2, order=TRUE)
    } else if (length(input$gene == 1)) {
      FeaturePlot(tj, features = c(input$gene), reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
    }
  }
  
  createSummary <- function(average) {
    tj$ovlp <- findOvlPCells(input$gene)
    Idents(tj) <- tj$ovlp
    
    tj$sample.ovlp <- paste(Idents(tj), tj$sample, sep = "_")
    Idents(tj) <- tj$sample.ovlp
    tj$sample.ovlp.cluster <- paste(Idents(tj), tj$seurat_clusters, sep = "_")
    
    Idents(tj) <- tj$sample
    tj$sample.cluster <- paste(Idents(tj), tj$seurat_clusters, sep = "_")
    
    Idents(tj) <- tj$sample.ovlp
    num_tj_ovlp <- length(WhichCells(tj, idents = "overlapping_cells_TJ"))
    
    num_tj <- length(tj$sample[which(tj$sample == "TJ")])
    
    pct_tj_ovlp <- (num_tj_ovlp/num_tj) * 100
    
    cat(paste("Total Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_tj_ovlp, "\n", sep=""))
    cat(paste("% Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_tj_ovlp, 2), nsmall = 2), sep=""))
    
    cat(paste("----------------------------------", "\n"))
    Idents(object = tj) <- tj$seurat_clusters
    num_clusters <- as.numeric(tail(levels(tj@meta.data$seurat_clusters), n=1))
    for (i in 0:num_clusters) {
      num_tj_cells_cluster <- 0
      
      Idents(tj) <- tj$sample.ovlp.cluster
      try(num_tj_cells_cluster <- length(WhichCells(tj, idents = paste("overlapping_cells_TJ", i, sep="_"))), silent=TRUE)
      
      if (average == TRUE) {
        Idents(tj) <- tj$sample.cluster
        num_tj_cells_cluster <- (num_tj_cells_cluster/length(WhichCells(tj, idents = paste("TJ", i, sep="_")))) * 100
        
        cat(paste("% TJ Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", format(round(num_tj_cells_cluster, 2), nsmall = 2), "% \n\n", sep=""))
      } else {
        Idents(tj) <- tj$sample.cluster
        cat(paste("TJ Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_tj_cells_cluster, " (out of ", length(WhichCells(tj, idents = paste("TJ", i, sep="_"))), ") \n\n", sep=""))
      }
      
    }
  }
  
  # Server Functions
  output$umap_plot <- renderPlot({
    createOvlpPlot()
  })
  
  output$barplot <- renderPlot({
    createBarPlot(FALSE)
  })
  
  output$barplotAvg <- renderPlot({
    createBarPlot(TRUE)
  })
  
  output$summary <- renderPrint({
    if (input$gene != "") {
      createSummary(FALSE)
    }
  })
  
  output$summaryAvg <- renderPrint({
    if (input$gene != "") {
      createSummary(TRUE)
    }
  })
  
}
shinyApp(ui = ui, server = server)