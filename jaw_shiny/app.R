###########
# Jaw APP #
###########

# Install Packages if not installed
if (!requireNamespace("shiny", quietly = TRUE))
  install.packages("shiny")
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")
if (!requireNamespace("Matrix", quietly = TRUE))
  install.packages("Matrix")
if (!requireNamespace("reticulate", quietly = TRUE))
  install.packages("reticulate")
if (!requireNamespace("stringr", quietly = TRUE))
  install.packages("stringr")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

# Load Packages
library("shiny")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("ggplot2")

# Load Data
combined <- readRDS("data/jaw.rds")
gene_names <- rownames(combined@assays$RNA)
# jaw <- readRDS("C:/Users/miles/Downloads/d_tooth/tooth_scripts/jaw_shiny/data/jaw.rds")

# UI
ui <- fluidPage(
  
  # App title ----
  titlePanel("Paint Expression of Cells"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      selectizeInput(inputId = "gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 10))
      # textInput(inputId = "cluster", label = "Cluster to Split", value = "0", width = NULL),
      # textInput(inputId = "resolution", label = "Resolution", value = "0.2", width = NULL)
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(id = "tabs", type = "tabs",
        tabPanel("LJ vs UJ Overlap", value="bc_ovlp_plot", plotOutput("bc_ovlp_plot", width = "100%", height="500px")),
        tabPanel("LJ vs UJ", value="bcplot", plotOutput("bcplot", width = "100%", height="500px")),
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
      expr <- FetchData(object = combined, vars = gene)
      pos_cells <- c(colnames(combined[, which(x = expr > 1)]), pos_cells)
    }
    
    # Find the overlap of the positive cells
    num_genes <- length(genes)
    counts <- table(pos_cells)
    ovlp_cells <- rownames(counts[counts >= num_genes])
    combined <- SetIdent(combined, cells=ovlp_cells, value="overlapping_cells")
    combined <- SetIdent(combined, cells=setdiff(WhichCells(combined), ovlp_cells), value="non_overlapping_cells")
    combined$ovlp <- combined@active.ident
    return(combined$ovlp)
  }
  
  createBarPlot <- function(average) {
    combined$ovlp <- findOvlPCells(input$gene)
    Idents(combined) <- combined$ovlp
    test <- combined[, WhichCells(combined, idents = "overlapping_cells")]
    num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
    total_lj <- c()
    total_uj <- c()
    Idents(test) <- test$sample
    Idents(combined) <- combined$sample
    test$cond.cluster <- paste(Idents(test), test$seurat_clusters, sep = "_")
    combined$cond.cluster <- paste(Idents(combined), combined$seurat_clusters, sep = "_")
    Idents(combined) <- combined$cond.cluster
    Idents(test) <- test$cond.cluster
    for (i in 0:num_clusters) {
      lj_cells_in_cluster <- 0
      uj_cells_in_cluster <- 0
      
      try(lj_cells_in_cluster <- length(WhichCells(test, idents = paste("LJ", i, sep="_"))), silent=TRUE)
      try(uj_cells_in_cluster <- length(WhichCells(test, idents = paste("UJ", i, sep="_"))), silent=TRUE)
      
      all_lj_cells_in_cluster <- length(WhichCells(combined, idents = paste("LJ", i, sep="_")))
      all_uj_cells_in_cluster <- length(WhichCells(combined, idents = paste("UJ", i, sep="_")))
      
      if (average == TRUE) {
        total_lj <- c(total_lj, round( (lj_cells_in_cluster/all_lj_cells_in_cluster) * 100, 2))
        total_uj <- c(total_uj, round( (uj_cells_in_cluster/all_uj_cells_in_cluster) * 100, 2))
      } else {
        total_lj <- c(total_lj, lj_cells_in_cluster)
        total_uj <- c(total_uj, uj_cells_in_cluster)
      }
      
    }
    df <- data.frame(condition <- c(rep("LJ", length(total_lj)), rep("UJ", length(total_uj))),
                     cluster_num <- c(0:num_clusters, 0:num_clusters),
                     value <- c(total_lj, total_uj))
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
  
  createSummary <- function(average) {
    combined$ovlp <- findOvlPCells(input$gene)
    Idents(combined) <- combined$ovlp
    
    combined$sample.ovlp <- paste(Idents(combined), combined$sample, sep = "_")
    Idents(combined) <- combined$sample.ovlp
    combined$sample.ovlp.cluster <- paste(Idents(combined), combined$seurat_clusters, sep = "_")
    
    Idents(combined) <- combined$sample
    combined$sample.cluster <- paste(Idents(combined), combined$seurat_clusters, sep = "_")
    
    Idents(combined) <- combined$sample.ovlp
    num_LJ_ovlp <- length(WhichCells(combined, idents = "overlapping_cells_LJ"))
    num_UJ_ovlp <- length(WhichCells(combined, idents = "overlapping_cells_UJ"))
    
    num_LJ <- length(combined$sample[which(combined$sample == "LJ")])
    num_UJ <- length(combined$sample[which(combined$sample == "UJ")])
    
    pct_LJ_ovlp <- (num_LJ_ovlp/num_LJ) * 100
    pct_UJ_ovlp <- (num_UJ_ovlp/num_UJ) * 100
    
    cat(paste("Total Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_LJ_ovlp+num_UJ_ovlp, "\n", sep=""))
    cat(paste("LJ vs UJ Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_LJ_ovlp, " vs ", num_UJ_ovlp, "\n", sep=""))
    cat(paste("% LJ vs UJ Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_LJ_ovlp, 2), nsmall = 2), "% vs ", format(round(pct_UJ_ovlp, 2), nsmall = 2), "% \n\n", sep=""))
    
    cat(paste("LJ Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_LJ_ovlp, "\n", sep=""))
    cat(paste("UJ Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", num_UJ_ovlp, "\n\n", sep=""))
    
    cat(paste("% LJ Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_LJ_ovlp, 2), nsmall = 2), "% \n", sep=""))
    cat(paste("% UJ Cells Expressing ", paste(input$gene, collapse = ' and '), ": ", format(round(pct_UJ_ovlp, 2), nsmall = 2), "% \n\n", sep=""))
    
    cat(paste("----------------------------------", "\n"))
    Idents(object = combined) <- combined$seurat_clusters
    num_clusters <- as.numeric(tail(levels(combined@meta.data$seurat_clusters), n=1))
    for (i in 0:num_clusters) {
      num_LJ_cells_cluster <- 0
      num_UJ_cells_cluster <- 0
      
      Idents(combined) <- combined$sample.ovlp.cluster
      try(num_LJ_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_LJ", i, sep="_"))), silent=TRUE)
      try(num_UJ_cells_cluster <- length(WhichCells(combined, idents = paste("overlapping_cells_UJ", i, sep="_"))), silent=TRUE)
      
      if (average == TRUE) {
        Idents(combined) <- combined$sample.cluster
        num_LJ_cells_cluster <- (num_LJ_cells_cluster/length(WhichCells(combined, idents = paste("LJ", i, sep="_")))) * 100
        num_UJ_cells_cluster <- (num_UJ_cells_cluster/length(WhichCells(combined, idents = paste("UJ", i, sep="_")))) * 100
        
        cat(paste("% LJ Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", format(round(num_LJ_cells_cluster, 2), nsmall = 2), "% \n", sep=""))
        cat(paste("% UJ Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", format(round(num_UJ_cells_cluster, 2), nsmall = 2), "% \n\n", sep=""))
      } else {
        Idents(combined) <- combined$sample.cluster
        cat(paste("LJ Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_LJ_cells_cluster, " (out of ", length(WhichCells(combined, idents = paste("LJ", i, sep="_"))), ") \n", sep=""))
        cat(paste("UJ Cells Expressing ", paste(input$gene, collapse = ' and '), " in cluster ", i, ": ", num_UJ_cells_cluster, " (out of ", length(WhichCells(combined, idents = paste("UJ", i, sep="_"))), ") \n\n", sep=""))
      }
      
    }
  }
  
  createOvlpPlot <- function() {
    if (length(input$gene) > 1) {
      combined$ovlp <- findOvlPCells(input$gene)
      Idents(combined) <- combined$ovlp
      DimPlot(combined, reduction="umap", group.by = "ident", pt.size=2, order=TRUE)
    } else if (length(input$gene == 1)) {
      FeaturePlot(combined, features = c(input$gene), reduction = "umap", split.by = "cond", pt.size = 2, label=TRUE, order = TRUE)
    }
  }
  
  createPlot <- function(split) {
    combined@active.assay <- "RNA"
    if (split == "cond") {
      FeaturePlot(combined, features = c(input$gene), split.by = "cond", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
    } else {
      Idents(combined) <- combined$sample
      only_b1_c1 <- combined[,WhichCells(combined,idents = c("b1", "c1"))]
      Idents(only_b1_c1) <- only_b1_c1$seurat_clusters
      FeaturePlot(only_b1_c1, features = c(input$gene), split.by = "sample", reduction = "umap", pt.size = 2, label=TRUE, order = TRUE)
    }
  }
  
  # Server Functions
  output$bcplot <- renderPlot({
    createPlot("cond")
  })
  
  output$bc_ovlp_plot <- renderPlot({
    createOvlpPlot()
  })
  
  output$bcplot <- renderPlot({
    createPlot("cond")
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