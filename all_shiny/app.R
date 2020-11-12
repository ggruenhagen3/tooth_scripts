###########
# All APP #
###########
library("shiny")
library("Seurat")
library("Matrix")
library("reticulate")
library("stringr")
library("ggplot2", warn.conflicts = FALSE, verbose = FALSE, quietly = TRUE)
library("cowplot", warn.conflicts = FALSE, verbose = FALSE, quietly = TRUE)

# Save space
obj <- list()
obj[["mes"]] <- readRDS("data/mes.rds")
# obj[["mes"]]@assays$RNA@scale.data <- matrix()
# empty_index <- which(rowSums(as.matrix(obj[["mes"]]@assays$RNA@counts)) == 0)
# obj[["mes"]]@assays$RNA@counts <- obj[["mes"]]@assays$RNA@counts[-empty_index,]
# obj[["mes"]]@assays$RNA@data <- obj[["mes"]]@assays$RNA@data[-empty_index,]
# obj[["mes"]]@assays$RNA@meta.features <- data.frame()
# obj[["mes"]]@assays$RNA@counts <- matrix()

obj[["epi"]] <- readRDS("data/epi.rds")
# obj[["epi"]]@assays$RNA@scale.data <- matrix()
# empty_index <- which(rowSums(as.matrix(obj[["epi"]]@assays$RNA@counts)) == 0)
# obj[["epi"]]@assays$RNA@counts <- obj[["epi"]]@assays$RNA@counts[-empty_index,]
# obj[["epi"]]@assays$RNA@data <- obj[["epi"]]@assays$RNA@data[-empty_index,]
# obj[["epi"]]@assays$RNA@meta.features <- data.frame()
# obj[["epi"]]@assays$RNA@counts <- matrix()

obj[["jaw"]] <- readRDS("data/jaw.rds")
# obj[["jaw"]]@assays$RNA@scale.data <- matrix()
# empty_index <- which(rowSums(as.matrix(obj[["jaw"]]@assays$RNA@counts)) == 0)
# obj[["jaw"]]@assays$RNA@counts <- obj[["jaw"]]@assays$RNA@counts[-empty_index,]
# obj[["jaw"]]@assays$RNA@data <- obj[["jaw"]]@assays$RNA@data[-empty_index,]
# obj[["jaw"]]@assays$RNA@counts <- matrix()

obj[["tj"]]  <- readRDS("data/tj.rds")
# obj[["tj"]]@assays$RNA@scale.data <- matrix()
# empty_index <- which(rowSums(as.matrix(obj[["tj"]]@assays$RNA@counts)) == 0)
# obj[["tj"]]@assays$RNA@counts <- obj[["tj"]]@assays$RNA@counts[-empty_index,]
# obj[["tj"]]@assays$RNA@data <- obj[["tj"]]@assays$RNA@data[-empty_index,]
# obj[["tj"]]@assays$RNA@counts <- matrix()

# saveRDS(obj[["mes"]], file = "C:/Users/miles/Downloads/d_tooth/tooth_scripts/all_shiny/data/mes.rds")
# saveRDS(obj[["epi"]], file = "C:/Users/miles/Downloads/d_tooth/tooth_scripts/all_shiny/data/epi.rds")
# saveRDS(obj[["jaw"]], file = "C:/Users/miles/Downloads/d_tooth/tooth_scripts/all_shiny/data/jaw.rds")
# saveRDS(obj[["tj"]], file = "C:/Users/miles/Downloads/d_tooth/tooth_scripts/all_shiny/data/tj.rds")

# mzebra_mouse <- read.table("data/mzebra_mouse.txt", sep = "\t", header = TRUE)
all_mouse <- read.table("data/all_mouse.txt", sep = "\t", header = TRUE)
obj[["epi"]]$cond[is.na(obj[["epi"]]$cond)] <- "INJR"
cichlid_gene_names <- unique(c(rownames(obj[["jaw"]]@assays$RNA), rownames(obj[["tj"]]@assays$RNA)))
mes_gene_names <- rownames(obj[["mes"]]@assays$RNA)
epi_gene_names <- rownames(obj[["epi"]]@assays$RNA)
mouse_gene_names <- unique(c(mes_gene_names, epi_gene_names)) # this should be the union between mes and epi

ui <- fillPage(
  navbarPage("All Tooth Datasets Painting App",
             id = "navbar",
             tabPanel("All",
                      
                      # Sidebar layout with input and output definitions ----
                      sidebarLayout(
                        
                        # Sidebar panel for inputs ----
                        sidebarPanel(
                          textOutput(outputId = "closest_cichlid_gene"),
                          selectizeInput(inputId = "all_gene_mouse", label = "Mouse Gene Name", choices = NULL, selected = "Celsr1", multiple = TRUE, options = list(maxOptions = 1000)),
                          textOutput(outputId = "closest_mouse_gene"),
                          selectizeInput(inputId = "all_gene_cichlid", label = "Cichlid Gene Name", choices = NULL, selected = "celsr1a", multiple = TRUE, options = list(maxOptions = 1000)),
                          width = 3
                        ),
                        
                        # Main panel for displaying outputs ----
                        mainPanel(
                          tags$style(type="text/css", '#mainPanel { height = 100%;}'),
                          id = "mainPanel",
                          plotOutput("all_umap_plot", width = "100%", height="500px"),
                          width = 9
                        )
                      )
                    ),
             tabPanel("Mouse Mesenchyme",
                      sidebarLayout(
                        sidebarPanel(
                          selectizeInput(inputId = "mes_gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 1000)),
                          width = 3
                        ),
                        
                        mainPanel(
                          
                          tabsetPanel(id = "mes_tabs", type = "tabs",
                                      tabPanel("UMAP", value="umap_plot", plotOutput("mes_umap_plot", width = "100%", height="500px")),
                                      tabPanel("UMAP By Condition", value="umap_plot", plotOutput("mes_umap_cond_plot", width = "100%", height="500px")),
                                      tabPanel("BarPlot Summary", value="mes_barplot", plotOutput("mes_barplot", width = "100%", height="500px")),
                                      tabPanel("BarPlot %", value="mes_barplotAvg", plotOutput("mes_barplotAvg", width = "100%", height="500px")),
                                      tabPanel("Summary", value="mes_summary", verbatimTextOutput("mes_summary")),
                                      tabPanel("Summary %", value="mes_summaryAvg", verbatimTextOutput("mes_summaryAvg"))
                          ),
                          width = 9
                        )
                      )),
             tabPanel("Mouse Epithelium",
                      sidebarLayout(
                        sidebarPanel(
                          selectizeInput(inputId = "epi_gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 1000)),
                          width = 3
                        ),
                        
                        mainPanel(
                          
                          tabsetPanel(id = "epi_tabs", type = "tabs",
                                      tabPanel("UMAP", value="umap_plot", plotOutput("epi_umap_plot", width = "100%", height="500px")),
                                      tabPanel("UMAP By Condition", value="umap_plot", plotOutput("epi_umap_cond_plot", width = "100%", height="500px")),
                                      tabPanel("BarPlot Summary", value="epi_barplot", plotOutput("epi_barplot", width = "100%", height="500px")),
                                      tabPanel("BarPlot %", value="epi_barplotAvg", plotOutput("epi_barplotAvg", width = "100%", height="500px")),
                                      tabPanel("Summary", value="epi_summary", verbatimTextOutput("epi_summary")),
                                      tabPanel("Summary %", value="epi_summaryAvg", verbatimTextOutput("epi_summaryAvg"))
                          ),
                          width = 9
                        )
                      )),
             tabPanel("Cichlid Jaw",
                      sidebarLayout(
                        sidebarPanel(
                          selectizeInput(inputId = "jaw_gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 1000)),
                          width = 3
                        ),
                        
                        mainPanel(
                          
                          tabsetPanel(id = "jaw_tabs", type = "tabs",
                                      tabPanel("UMAP", value="umap_plot", plotOutput("jaw_umap_plot", width = "100%", height="500px")),
                                      tabPanel("BarPlot Summary", value="jaw_barplot", plotOutput("jaw_barplot", width = "100%", height="500px")),
                                      tabPanel("BarPlot %", value="jaw_barplotAvg", plotOutput("jaw_barplotAvg", width = "100%", height="500px")),
                                      tabPanel("Summary", value="jaw_summary", verbatimTextOutput("jaw_summary")),
                                      tabPanel("Summary %", value="jaw_summaryAvg", verbatimTextOutput("jaw_summaryAvg"))
                          ),
                          width = 9
                        )
                      )),
             tabPanel("Cichlid Tooth",
                      sidebarLayout(
                        sidebarPanel(
                          selectizeInput(inputId = "tj_gene", label = "Gene Name", choices = NULL, selected = "fosb", multiple = TRUE, options = list(maxOptions = 1000)),
                          width = 3
                        ),
                        
                        mainPanel(
                          
                          tabsetPanel(id = "tj_tabs", type = "tabs",
                                      tabPanel("UMAP", value="tj_umap_plot", plotOutput("tj_umap_plot", width = "100%", height="500px")),
                                      tabPanel("BarPlot Summary", value="tj_barplot", plotOutput("tj_barplot", width = "100%", height="500px")),
                                      tabPanel("BarPlot %", value="tj_barplotAvg", plotOutput("tj_barplotAvg", width = "100%", height="500px")),
                                      tabPanel("Summary", value="tj_summary", verbatimTextOutput("tj_summary")),
                                      tabPanel("Summary %", value="tj_summaryAvg", verbatimTextOutput("tj_summaryAvg"))
                          ),
                          width = 9
                        )
                      ))
  )
)

server <- function(input, output, session) {
  # Helpful Input Suggestions
  updateSelectizeInput(session, "all_gene_mouse",   choices = mouse_gene_names, server = TRUE)
  updateSelectizeInput(session, "all_gene_cichlid", choices = cichlid_gene_names, server = TRUE)
  updateSelectizeInput(session, "mes_gene", choices = mes_gene_names, server = TRUE)
  updateSelectizeInput(session, "epi_gene", choices = epi_gene_names, server = TRUE)
  updateSelectizeInput(session, "jaw_gene", choices = cichlid_gene_names, server = TRUE)
  updateSelectizeInput(session, "tj_gene",  choices = cichlid_gene_names, server = TRUE)
  
  # Helper Functions
  findOvlPCells <- function(obj_str, genes) {
    clean_genes <- c()
    pos_cells <- c()
    for (gene in genes) {
      clean_genes <- genes
      expr <- FetchData(object = obj[[obj_str]], vars = gene)
      pos_cells <- c(colnames(obj[[obj_str]][, which(x = expr > 0)]), pos_cells)
    }
    
    # Find the overlap of the positive cells
    num_genes <- length(genes)
    counts <- table(pos_cells)
    ovlp_cells <- rownames(counts)[which(counts >= num_genes)]
    obj[[obj_str]] <- SetIdent(obj[[obj_str]], cells=ovlp_cells, value="overlapping_cells")
    obj[[obj_str]] <- SetIdent(obj[[obj_str]], cells=setdiff(WhichCells(obj[[obj_str]]), ovlp_cells), value="non_overlapping_cells")
    obj[[obj_str]]$ovlp <- obj[[obj_str]]@active.ident
    return(obj[[obj_str]]$ovlp)
  }
  
  createBarPlot <- function(obj_str, average, genes) {
    obj[[obj_str]]$ovlp <- findOvlPCells(obj_str, genes)
    Idents(obj[[obj_str]]) <- obj[[obj_str]]$ovlp
    conditions <- unique(obj[[obj_str]]$cond)
    ovlp_cells <- WhichCells(obj[[obj_str]], idents = "overlapping_cells")
    df <- data.frame()
    num_clusters <- as.numeric(tail(levels(obj[[obj_str]]@meta.data$seurat_clusters), n=1))
    for (i in 0:num_clusters) {
      Idents(obj[[obj_str]]) <- obj[[obj_str]]$seurat_clusters
      this_cluster <- WhichCells(obj[[obj_str]], idents = i)
      
      Idents(obj[[obj_str]]) <- obj[[obj_str]]$cond
      for (cond in conditions) {
        this_cond <- WhichCells(obj[[obj_str]], idents = cond)
        stat_to_display <- length(ovlp_cells[which(ovlp_cells %in% this_cluster & ovlp_cells %in% this_cond)])
        if(average == TRUE) {
          stat_to_display <- round( (stat_to_display / length(this_cond[which(this_cond %in% this_cluster)]))*100 ,2)
        }
        newRow <- t(c(cond, i, stat_to_display))
        df <- rbind(df, newRow)
      } # end cond for
    } # end cluster for
    colnames(df) <- c("condition", "cluster_num", "value")
    df$value <- as.numeric(as.vector(df$value))
    df$condition <- as.character(df$condition)
    df$cluster_num <- factor(df$cluster_num, levels = 0:num_clusters)
    my_title <- paste("Number of Cells Expressing", paste(genes, collapse = ' and '), "per Cluster")
    my_ylab <- "Number of Cells"
    if (average == TRUE) {
      my_title <- paste("% Cells Expressing", paste(genes, collapse = ' and '), "per Cluster")
      my_ylab <- "% Cells"
    }
    p <- ggplot(df, aes(fill=condition, x=cluster_num, y=value)) +
      geom_bar(position="dodge", stat="identity") +
      theme_minimal() +
      ggtitle(my_title) +
      xlab("Cluster") +
      ylab(my_ylab) +
      geom_text(aes(label=value), vjust=1.6, color="black", position = position_dodge(0.9), size=3.5)
    theme_minimal()
    p
  }
  
  createOvlpPlot <- function(obj_str, genes, split) {
    if (length(genes) > 1) {
      obj[[obj_str]]$ovlp <- findOvlPCells(obj_str, genes)
      Idents(obj[[obj_str]]) <- obj[[obj_str]]$ovlp
      DimPlot(obj[[obj_str]], reduction="umap", split.by = split, group.by = "ident", pt.size=2, order=TRUE)
    } else if (length(genes == 1)) {
      FeaturePlot(obj[[obj_str]], features = c(genes), reduction = "umap", split.by = split, pt.size = 2, label=TRUE, order = TRUE)
    }
  }
  
  createSummary <- function(obj_str, average, genes) {
    obj[[obj_str]]$ovlp <- findOvlPCells(obj_str, genes)
    Idents(obj[[obj_str]]) <- obj[[obj_str]]$ovlp
    ovlp_cells <- WhichCells(obj[[obj_str]], idents = "overlapping_cells")

    num_obj_ovlp  <- length(ovlp_cells)
    num_all_cells <- ncol(obj[[obj_str]])
    
    pct_obj_ovlp <- (num_obj_ovlp/num_all_cells) * 100
    
    cat(paste("Total Cells Expressing ", paste(genes, collapse = ' and '), ": ", num_obj_ovlp, "\n", sep=""))
    cat(paste("% Cells Expressing ", paste(genes, collapse = ' and '), ": ", format(round(pct_obj_ovlp, 2), nsmall = 2), "\n", sep=""))
    
    cat(paste("----------------------------------", "\n"))
    Idents(object = obj[[obj_str]]) <- obj[[obj_str]]$seurat_clusters
    num_clusters <- as.numeric(tail(levels(obj[[obj_str]]@meta.data$seurat_clusters), n=1))
    for (i in 0:num_clusters) {
      this_cells <- WhichCells(obj[[obj_str]], idents = i)
      num_obj_cells_cluster <- length(ovlp_cells[which(ovlp_cells %in% this_cells)])
      
      if (average == TRUE) {
        num_obj_cells_cluster <- (num_obj_cells_cluster/length(this_cells)) * 100
        cat(paste("% Cells Expressing ", paste(genes, collapse = ' and '), " in cluster ", i, ": ", format(round(num_obj_cells_cluster, 2), nsmall = 2), "% \n\n", sep=""))
      } else {
        cat(paste("Cells Expressing ", paste(genes, collapse = ' and '), " in cluster ", i, ": ", num_obj_cells_cluster, " (out of ", length(this_cells), ") \n\n", sep=""))
      }
      
    }
  }
  
  createAllUmapPlot <- function(mouse_gene, cichlid_gene) {
    mouse_gene <- as.character(mouse_gene)
    cichlid_gene <- as.character(cichlid_gene)
    myplots <- list()
    if (length(mouse_gene) > 0 && mouse_gene %in% rownames(obj[["mes"]])) {
      myplots[[length(myplots)+1]] <- FeaturePlot(obj[["mes"]], mouse_gene, reduction = "umap", order = TRUE, pt.size=1.5, label = TRUE)
    }
    if (length(mouse_gene) >0 && mouse_gene %in% rownames(obj[["epi"]])) {
      myplots[[length(myplots)+1]] <- FeaturePlot(obj[["epi"]], mouse_gene, reduction = "umap", order = TRUE, pt.size=1.5, label = TRUE)
    }
    if (length(cichlid_gene) > 0 && cichlid_gene %in% rownames(obj[["jaw"]])) {
      myplots[[length(myplots)+1]] <- FeaturePlot(obj[["jaw"]], cichlid_gene, reduction = "umap", order = TRUE, pt.size=1.5, label = TRUE)
    }
    if (length(cichlid_gene) > 0 && cichlid_gene %in% rownames(obj[["tj"]])) {
      myplots[[length(myplots)+1]] <- FeaturePlot(obj[["tj"]],  cichlid_gene, reduction = "umap", order = TRUE, pt.size=1.5, label = TRUE)
    }
    plot_grid(plotlist= myplots)
  }
  
  findClosestCichlid <- function(mouse_gene) {
    this_cichlid_genes <- all_mouse[which(all_mouse[,2] == mouse_gene),1]
    if (length(this_cichlid_genes) > 1) {
      upper_cichlid_gene <- this_cichlid_genes[which(startsWith(tolower(this_cichlid_genes), tolower(mouse_gene)))]
      if (length(upper_cichlid_gene) == 1) {
        cichlid_gene <- upper_cichlid_gene
      } else {
        cichlid_gene <- this_cichlid_genes[1]
      } # end bad multiple
    } else {
      cichlid_gene <- this_cichlid_genes
    }
    
    return(cichlid_gene)
  }
  
  findClosestMouse <- function(cichlid_gene) {
    this_mouse_genes <- all_mouse[which(all_mouse[,1] == cichlid_gene),2]
    if (length(this_mouse_genes) > 1) {
      upper_mouse_gene <- this_rows[which(startsWith(tolower(this_mouse_genes), tolower(cichlid_gene)))]
      if (length(upper_mouse_gene) == 1) {
        mouse_gene <- upper_mouse_gene
      } else {
        mouse_gene <- this_mouse_genes[1]
      } # end bad multiple
    } else {
      mouse_gene <- this_mouse_genes
    }
    
    return(mouse_gene)
  }
  
  # Server Functions
  output$all_umap_plot <- renderPlot({
    if (length(input$all_gene_cichlid) > 0 && length(input$all_gene_mouse) == 0) {
      mouse_gene <- findClosestMouse(input$all_gene_cichlid)
      createAllUmapPlot(mouse_gene, input$all_gene_cichlid)
    } else if (length(input$all_gene_cichlid) == 0 && length(input$all_gene_mouse) > 0) {
      cichlid_gene <- findClosestCichlid(input$all_gene_mouse)
      createAllUmapPlot(input$all_gene_mouse, cichlid_gene)
    } else if (length(input$all_gene_cichlid) == 1 && length(input$all_gene_mouse) == 1) {
      createAllUmapPlot(input$all_gene_mouse, input$all_gene_cichlid)
    }
  })
  
  output$closest_cichlid_gene <- renderText({
    if (length(input$all_gene_mouse) > 0) {
      str <- "Closest Cichlid Gene: "
      str <- paste(str, as.character(findClosestCichlid(input$all_gene_mouse)))
      return(str)
    }
  })
  
  output$closest_mouse_gene <- renderText({
    if (length(input$all_gene_cichlid) > 0) {
      str <- "Closest Mouse Gene: "
      str <- paste(str, as.character(findClosestMouse(input$all_gene_cichlid)))
      return(str)
    }
  })
  
  output$mes_umap_plot <- renderPlot({
    if (length(input$mes_gene) != 0) {
      genes <- input$mes_gene
      createOvlpPlot("mes", genes, NULL)
    }
  })
  
  output$mes_umap_cond_plot <- renderPlot({
    if (length(input$mes_gene) != 0) {
      genes <- input$mes_gene
      createOvlpPlot("mes", genes, "cond")
    }
  })
  
  output$mes_barplot <- renderPlot({
    if (length(input$mes_gene) != 0) {
      genes <- input$mes_gene
      createBarPlot("mes", FALSE, genes)
    }
  })
  
  output$mes_barplotAvg <- renderPlot({
    if (length(input$mes_gene) != 0) {
      genes <- input$mes_gene
      createBarPlot("mes", TRUE, genes)
    }
  })
  
  output$mes_summary <- renderPrint({
    if (length(input$mes_gene) != 0) {
      genes <- input$mes_gene
      createSummary("mes", FALSE, genes)
    }
  })
  
  output$mes_summaryAvg <- renderPrint({
    if (length(input$mes_gene) != 0) {
      genes <- input$mes_gene
      createSummary("mes", TRUE, genes)
    }
  })
  
  output$epi_umap_plot <- renderPlot({
    if (length(input$epi_gene) != 0) {
      genes <- input$epi_gene
      createOvlpPlot("epi", genes, NULL)
    }
  })
  
  output$epi_umap_cond_plot <- renderPlot({
    if (length(input$epi_gene) != 0) {
      genes <- input$epi_gene
      createOvlpPlot("epi", genes, "cond")
    }
  })
  
  output$epi_barplot <- renderPlot({
    if (length(input$epi_gene) != 0) {
      genes <- input$epi_gene
      createBarPlot("epi", FALSE, genes)
    }
  })
  
  output$epi_barplotAvg <- renderPlot({
    if (length(input$epi_gene) != 0) {
      genes <- input$epi_gene
      createBarPlot("epi", TRUE, genes)
    }
  })
  
  output$epi_summary <- renderPrint({
    if (length(input$epi_gene) != 0) {
      genes <- input$epi_gene
      createSummary("epi", FALSE, genes)
    }
  })
  
  output$epi_summaryAvg <- renderPrint({
    if (length(input$epi_gene) != 0) {
      genes <- input$epi_gene
      createSummary("epi", TRUE, genes)
    }
  })
  
  output$jaw_umap_plot <- renderPlot({
    if (length(input$jaw_gene) != 0) {
      genes <- input$jaw_gene
      createOvlpPlot("jaw", genes, NULL)
    }
  })
  
  output$jaw_barplot <- renderPlot({
    if (length(input$jaw_gene) != 0) {
      genes <- input$jaw_gene
      createBarPlot("jaw", FALSE, genes)
    }
  })
  
  output$jaw_barplotAvg <- renderPlot({
    if (length(input$jaw_gene) != 0) {
      genes <- input$jaw_gene
      createBarPlot("jaw", TRUE, genes)
    }
  })
  
  output$jaw_summary <- renderPrint({
    if (length(input$jaw_gene) != 0) {
      genes <- input$jaw_gene
      createSummary("jaw", FALSE, genes)
    }
  })
  
  output$jaw_summaryAvg <- renderPrint({
    if (length(input$jaw_gene) != 0) {
      genes <- input$jaw_gene
      createSummary("jaw", TRUE, genes)
    }
  })
  
  output$tj_umap_plot <- renderPlot({
    if (length(input$tj_gene) != 0) {
      genes <- input$tj_gene
      createOvlpPlot(obj[["tj"]], genes, NULL)
    }
  })
  
  output$tj_barplot <- renderPlot({
    if (length(input$tj_gene) != 0) {
      genes <- input$tj_gene
      createBarPlot(obj[["tj"]], FALSE, genes)
    }
  })
  
  output$tj_barplotAvg <- renderPlot({
    if (length(input$tj_gene) != 0) {
      genes <- input$tj_gene
      createBarPlot(obj[["tj"]], TRUE, genes)
    }
  })
  
  output$tj_summary <- renderPrint({
    if (length(input$tj_gene) != 0) {
      genes <- input$tj_gene
      createSummary(obj[["tj"]], FALSE, genes)
    }
  })
  
  output$tj_summaryAvg <- renderPrint({
    if (length(input$tj_gene) != 0) {
      genes <- input$tj_gene
      createSummary(obj[["tj"]], TRUE, genes)
    }
  })

  
}
shinyApp(ui = ui, server = server)