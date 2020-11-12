#############
# JPOOL APP #
#############
library("shiny")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("ggplot2")
jpool <- readRDS("data/jpool.rds")
gene_names <- rownames(jpool@assays$RNA)

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
      expr <- FetchData(object = jpool, vars = gene)
      pos_cells <- c(colnames(jpool[, which(x = expr > 1)]), pos_cells)
    }
    
    # Find the overlap of the positive cells
    num_genes <- length(genes)
    counts <- table(pos_cells)
    ovlp_cells <- rownames(counts[counts >= num_genes])
    jpool <- SetIdent(jpool, cells=ovlp_cells, value="overlapping_cells")
    jpool <- SetIdent(jpool, cells=setdiff(WhichCells(jpool), ovlp_cells), value="non_overlapping_cells")
    jpool$ovlp <- jpool@active.ident
    return(jpool$ovlp)
  }
  
  createBarPlot <- function(average) {
    jpool$ovlp <- findOvlPCells(input$gene)
    Idents(jpool) <- jpool$ovlp
    # jpool$sample <- rep("JPOOL", ncol(jpool))
    test <- jpool[, WhichCells(jpool, idents = "overlapping_cells")]
    num_clusters <- as.numeric(tail(levels(jpool@meta.data$seurat_clusters), n=1))
    total_jpool <- c()
    Idents(test) <- test$sample
    Idents(jpool) <- jpool$sample
    test$cond.cluster <- paste(Idents(test), test$seurat_clusters, sep = "_")
    jpool$cond.cluster <- paste(Idents(jpool), jpool$seurat_clusters, sep = "_")
    Idents(jpool) <- jpool$cond.cluster
    Idents(test) <- test$cond.cluster
    for (i in 0:num_clusters) {
      jpool_cells_in_cluster <- 0
      try(jpool_cells_in_cluster <- length(WhichCells(test, idents = paste("JPOOL", i, sep="_"))), silent=TRUE)
      all_jpool_cells_in_cluster <- length(WhichCells(jpool, idents = paste("JPOOL", i, sep="_")))
      
      if (average == TRUE) {
        total_jpool <- c(total_jpool, round( (jpool_cells_in_cluster/all_jpool_cells_in_cluster) * 100, 2))
      } else {
        total_jpool <- c(total_jpool, jpool_cells_in_cluster)
      }
      
    }
    df <- data.frame(condition <- c(rep("JPOOL", length(total_jpool))),cluster_num <- c(0:num_clusters),
                     value <- c(total_jpool))
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
      jpool$ovlp <- findOvlPCells(input$gene)
      Idents(jpool) <- jpool$ovlp
      DimPlot(jpool, reduction="umap", group.by = "ident", pt.size=2, order=TRUE)
    } else if (length(input$gene == 1)) {
      FeaturePlot(jpool, features = c(input$gene), reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
    }
  }
  
  createSummary <- function(average) {
    jpool$ovlp <- findOvlPCells(input$gene)
    Idents(jpool) <- jpool$ovlp
    
    jpool$sample.ovlp <- paste(Idents(jpool), jpool$sample, sep = "_")
    Idents(jpool) <- jpool$sample.ovlp
    jpool$sample.ovlp.cluster <- paste(Idents(jpool), jpool$seurat_clusters, sep = "_")
    
    Idents(jpool) <- jpool$sample
    jpool$sample.cluster <- paste(Idents(jpool), jpool$seurat_clusters, sep = "_")
    
    Idents(jpool) <- jpool$sample.ovlp
    num_jpool_ovlp <- length(WhichCells(jpool, idents = "overlapping_cells_JPOOL"))
    
    num_jpool <- length(jpool$sample[which(jpool$sample == "JPOOL")])
    
    pct_jpool_ovlp <- (num_jpool_ovlp/num_jpool) * 100
    
    cat(paste("Total Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_jpool_ovlp, "\n", sep=""))
    cat(paste("% Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_jpool_ovlp, 2), nsmall = 2), sep=""))
    
    cat(paste("----------------------------------", "\n"))
    Idents(object = jpool) <- jpool$seurat_clusters
    num_clusters <- as.numeric(tail(levels(jpool@meta.data$seurat_clusters), n=1))
    for (i in 0:num_clusters) {
      num_jpool_cells_cluster <- 0
      
      Idents(jpool) <- jpool$sample.ovlp.cluster
      try(num_jpool_cells_cluster <- length(WhichCells(jpool, idents = paste("overlapping_cells_JPOOL", i, sep="_"))), silent=TRUE)
      
      if (average == TRUE) {
        Idents(jpool) <- jpool$sample.cluster
        num_jpool_cells_cluster <- (num_jpool_cells_cluster/length(WhichCells(jpool, idents = paste("JPOOL", i, sep="_")))) * 100
        
        cat(paste("% jpool Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", format(round(num_jpool_cells_cluster, 2), nsmall = 2), "% \n\n", sep=""))
      } else {
        Idents(jpool) <- jpool$sample.cluster
        cat(paste("jpool Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_jpool_cells_cluster, " (out of ", length(WhichCells(jpool, idents = paste("JPOOL", i, sep="_"))), ") \n\n", sep=""))
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