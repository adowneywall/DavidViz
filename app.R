#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(tidyverse)
library(ggthemes)
library(circlize)
library(cowplot)
library(BiocManager)
#library(devtools)
options(repos = BiocManager::repositories())
#BiocManager::install('ComplexHeatmap')
# bioc_ls <- installed.packages() %>%
#   rownames() %>%
#   devtools::package_info() %>%
#   filter(source == "Bioconductor") %>%
#   select(package)
# if(!"ComplexHeatmap" %in% bioc_ls$package){
#   options(repos = BiocManager::repositories())
#   BiocManager::install('ComplexHeatmap')
# }
library(ComplexHeatmap)

#### Similarity Calculations ####

## Version 1: Calculates similarity as the proportion of genes shared relative to the smaller of the two gene sets
simProp_min <- function(x,y){
  z <- sum(x %in% y)/min(c(length(x),length(y)))
  return(z)
}

## Version 2: Calculates similarity as proportion of genes shared relative to the total number of unique genes across all gene sets
simProp_max <- function(x,y){
  m <- length(unique(c(x,y)))
  z <- sum(x %in% y)/m
  return(z)
}

jaccard <- function(x,y) {
  intersection = length(intersect(x,y))
  union = length(x) + length(y) - intersection
  return (intersection/union)
}

#### Similarity Matrix ####

# Generates pairwise similarity matrix when provided list of gene sets. Optional argument - matrix labels 
# (should be in the same order as the names of your gene sets).

simMatrix <- function(g_list,labels=NULL,sim_method='simProp_min'){
  gene_sim_mat <- matrix(ncol=length(g_list),nrow=length(g_list),data = 0)
  for(i in 1:length(g_list)){
    if(sim_method == 'simProp_max'){
      gene_sim_mat[,i] <- unlist(lapply(g_list,simProp_max,y=g_list[[i]]))
    }else{
      if(sim_method == 'jaccard'){
        gene_sim_mat[,i] <- unlist(lapply(g_list,jaccard,y=g_list[[i]]))
      }else{
        gene_sim_mat[,i] <- unlist(lapply(g_list,simProp_min,y=g_list[[i]]))
      }
    } 
  }
  colnames(gene_sim_mat) <- labels
  rownames(gene_sim_mat) <- labels
  return(gene_sim_mat)
}

#### Term Matrix #####
term_matrix <- function(x,viz_variable='FDR'){
  if(viz_variable == 'FDR'){
    # FDR
    out  <- x %>%  
      mutate(log_FDR= -log10(Top_FDR),
             lnc_term = paste0(GeneSetID,"_",Top_Term_Name)) %>%
      select(c('GeneSetID','Top_Term_Name','log_FDR','lnc_term')) %>%
      distinct(lnc_term,.keep_all=T) %>%
      select(c('GeneSetID','Top_Term_Name','log_FDR')) %>%
      pivot_wider(names_from = GeneSetID,values_from = log_FDR,values_fill = 0)
  }else{
    # Cluster Enrichment Score
    out <- x %>%  
      mutate(lnc_term = paste0(GeneSetID,"_",Top_Term_Name)) %>%
      select(c('GeneSetID','Top_Term_Name','Enrich_Score','lnc_term')) %>%
      distinct(lnc_term,.keep_all=T) %>%
      select(c('GeneSetID','Top_Term_Name','Enrich_Score')) %>%
      pivot_wider(names_from = GeneSetID,values_from = Enrich_Score,values_fill = 0)
  }
  
  out_lab <- out$Top_Term_Name
  out <- data.frame(out[,-1])
  out <- as.matrix(out)
  rownames(out) <- out_lab 
  return(out)
}

##### Default Heat Map #####
default_heatmap <- function(out,clusterMethod='cdef',
                            colorScale='Fixed',g_mat=NULL,legend_name='-log10(FDR)',
                            row_font=8,column_font=8){
  if(colorScale == 'Fixed'){
    if(legend_name == '-log10(FDR)'){
      col_fun = colorRamp2(c(0,0.1,20,40), c("white",'chartreuse2','coral','darkorchid2'))
    }else{
      col_fun = colorRamp2(c(0,0.1,10,20), c("white",'chartreuse2','coral','darkorchid2'))
    }
  }else{
    col_fun = colorRamp2(c(0,0.1,max(out)/2,max(out)+1), c("white",'chartreuse2','coral','darkorchid2'))
  }
  
  if(clusterMethod == 'cdef'){
    Heatmap(out,
            col = col_fun,
            rect_gp = gpar(col= "white"),
            #cluster_rows = function(f) hclust(dist(f)),
            clustering_distance_rows = function(x, y) 1 - cor(x, y),
            row_names_side = "left", row_dend_side = "right",
            row_dend_width = unit(5, "cm"),
            column_dend_height = unit(1, "cm"),
            row_names_gp = gpar(fontsize = row_font),
            column_names_gp = gpar(fontsize = column_font),
            row_names_max_width = max_text_width(
              rownames(out), 
              gp = gpar(fontsize = 12)
            ),
            heatmap_legend_param = list(
              title = legend_name))
  }else{
    Heatmap(out,
            col = col_fun,
            rect_gp = gpar(col= "white"),
            cluster_rows = hclust(dist(g_mat)),
            row_names_side = "left", row_dend_side = "right",
            row_dend_width = unit(5, "cm"),
            column_dend_height = unit(1, "cm"),
            row_names_gp = gpar(fontsize = row_font),
            column_names_gp = gpar(fontsize = column_font),
            row_names_max_width = max_text_width(
              rownames(out), 
              gp = gpar(fontsize = 12)
            ),
            heatmap_legend_param = list(
              title = legend_name))
  }
}

#### Geneset Heat Map ####
gene_set_heatmap <- function(x){
  col_fun = colorRamp2(c(0,0.5,0.9999,1), c('chartreuse2','coral','darkorchid2','white'))
  Heatmap(x,
          col = col_fun,
          rect_gp = gpar(col= "white"),
          #cluster_rows = function(f) hclust(dist(f)),
          #clustering_distance_rows = function(x, y) 1 - cor(x, y),
          row_names_side = "left", row_dend_side = "right",
          row_dend_width = unit(3, "cm"),
          column_dend_height = unit(3, "cm"),
          row_names_gp = gpar(fontsize = 9),
          row_names_max_width = max_text_width(
            rownames(x), 
            gp = gpar(fontsize = 10)
          ),
          column_names_gp = gpar(fontsize = 9),
          column_names_max_height = max_text_height(
            colnames(x), 
            gp = gpar(fontsize = 500)
          ),
          heatmap_legend_param = list(
            title = 'Dissimilarity'))
}

#### Enrichment Visualization Function #####
enrich_viz <- function(df_d=df_d,
                       group=group,
                       CLUSTER_SCORE =  3,
                       FDR_SCORE = 0.001,
                       GOTERM = "BP",
                       CLUSTER = 'cdef'){
  df_red <- df_d %>%
    filter(Score >= CLUSTER_SCORE,
           Top_FDR <= FDR_SCORE,
           grepl(paste0(GOTERM,collapse = "|"),Top_Category)) %>%
    arrange(desc(Score)) %>%
    mutate(Enrich_Score = Score) %>%
    relocate(Enrich_Score,.after = 'GeneSetID')
  
  group <- group %>%
    filter(grepl(paste0(GOTERM,collapse = "|"),Term_Type))
  
  ##### Generate term  matrices #####
  out_FDR <- term_matrix(df_red,viz_variable = 'FDR')
  #out_ES <- term_matrix(df_red,viz_variable = 'Enrichment_Score') # Removed this functionality for time being
  
  ### Generate gene set similarity matrix for each unique annotation cluster in file
  # Generates a heatmap to visualize these similarities.
  # Similarity based on the proportion of genes shared between each unique DAVID annotation
  # Keep only unique clusters
  df_red_unique <- df_red %>%
    distinct(Top_Term_Name,.keep_all = T)
  
  g_list <- str_split(df_red$All_Unique_Genes,pattern = ',')
  g_list_unique <- unique(df_red$Top_Term_Name)
  
  list_unique <- function(x,y,w){
    return(unique(unlist(w[which(y$Top_Term_Name == x)])))
  }
  out <- sapply(g_list_unique,list_unique,y=df_red,w=g_list)

  temp <- str_split(group$Term,pattern='~',simplify = T)
  group_red <- group[toupper(temp[,2]) %in% toupper(g_list_unique),]
  g_list_unique_DBoverlap <- g_list_unique[toupper(g_list_unique) %in% toupper(temp[toupper(temp[,2]) %in% toupper(g_list_unique),2])]
  out_list_DBoverlap <- out[toupper(g_list_unique) %in% toupper(temp[toupper(temp[,2]) %in% toupper(g_list_unique),2])]
  
  ### Genes present in genesets from DAVID file ###
  if(CLUSTER == 'cdef'){
    g_mat = NULL
  }else{
    if(CLUSTER == 'cmin'){
      g_mat = simMatrix(g_list = out,
                        labels = g_list_unique,
                        sim_method = 'simProp_min')
      # g_mat2 = simMatrix(g_list = group_red$Gene_list,
      #                    labels = toupper(temp[toupper(temp[,2]) %in% toupper(g_list_unique),2]))
    }else{
      if(CLUSTER == 'cmax'){
        g_mat = simMatrix(g_list = out,
                          labels = g_list_unique,
                          sim_method = 'simProp_min')
      }else{
        g_mat = simMatrix(g_list = out,
                          labels = g_list_unique,
                          sim_method = 'jaccard')
      }
    }
  }
  
  #gene_set_heatmap(1-g_mat_V1)
  return(list(out_FDR,g_mat))
}
  
#### Shiny Code ####
# Define UI for application that draws a histogram
ui <- fluidPage(
  useShinyjs(),
  # Application title
  titlePanel("DAVID Heatmap Analysis"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      width = 2,
      fileInput(inputId = "david_file",buttonLabel = "Upload",label="Upload enrichment file"),
      helpText("Uploaded tab-delimited file must be in the standard formated generated by the DAVID aggregator function"),
      tags$a(href="https://github.com/adowneywall/Davidaggregator", 
             "Link to DAVID aggregator program"),
      helpText(""),
      ## Buttons and parameters
      numericInput(inputId = "CS",
                   label = "Minimum Cluster Score:",
                   value = 2),
      verbatimTextOutput("CS_error"),
      numericInput(inputId = "FDR",
                   label = "FDR Threshold:",
                   value = 0.05,max = 1),
      verbatimTextOutput("FDR_error"),
      checkboxGroupInput(inputId = "Process",
                         label = "Go Terms:",
                         choices = c("Biological Processes (BP)" = "BP",
                                     "Molecular Functions (MF)" = "MF",
                                     "Cellular Components (CC)" = "CC")),
      verbatimTextOutput("process_error"),
      radioButtons(inputId = "cluster",
                   label = "Gene Ontology Clustering Method:",
                   choices = c("Enrichment metric (default)" = "cdef",
                               "Cluster Min" = "cmin",
                               "Cluster Max" = "cmax",
                               "Jaccard" = "cj")),
      helpText("Cluster Min : Calculates similarity as the proportion of genes shared relative to the smaller of the two gene sets"),
      helpText("Cluster Max : Calculates similarity as proportion of genes shared relative to the total number of unique genes across all gene sets"),
      sliderInput("height",
                  "Figure Height:",
                  min = 1,
                  max = 3000,
                  step = 10,
                  value = 1400),
      sliderInput("width",
                  "Figure Width:",
                  min = 1,
                  max = 3000,
                  step = 10,
                  value = 1400),
      sliderInput("font_r",
                  "Row Font:",
                  min = 1,
                  max = 30,
                  step = 0.1,
                  value = 12),
      sliderInput("font_c",
                  "Column Font:",
                  min = 1,
                  max = 30,
                  step = 0.1,
                  value = 12),
      
      actionButton("plot", "Plot"),
      actionButton("download", "Download"),
      verbatimTextOutput("error"),
      
      conditionalPanel(
        "false", # always hide the download button
        downloadButton(outputId = "downloadData", "Download tab-delimited (.tsv) file")
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("plot",width="100%")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  options(shiny.maxRequestSize=30*1024^2) 
  
  v <- reactiveValues(doPlot = FALSE,doDownload=FALSE)
  
  resulted_zip <- reactiveVal()
  main_dir <- getwd()
  
  # Reference File to read in
  ### All genes present in annotation on DAVID knowledgebase ###
  bp <- read_tsv('./ref/DAVID_BP_FAT_GOTERMS.tsv')
  mf <- read_tsv('./ref/DAVID_MF_FAT_GOTERMS.tsv')
  cc <- read_tsv('./ref/DAVID_CC_FAT_GOTERMS.tsv')
  up <- read_tsv('./ref/DAVID_UP_SEQ_GOTERMS.tsv')
  group <- rbind(bp,mf,cc,up)
  
  observeEvent(input$plot, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v$doPlot <- input$plot
    
    
    ## Generates plot in ui
    output$plot <- renderPlot({
      if (v$doPlot == FALSE) return()  
      
      ## Error Messages
      if(input$CS <= 0){
        output$CS_error <- renderText({
          validate(need(input$CS > 0, 'Minimum enrichment score needs to be greater than 0'))
        })
        return()
      }

      if(input$FDR <= 0){
        output$FDR_error <- renderText({
          validate(need(input$FDR > 0, 'FDR threshold needs to be greater than 0'))
        })
        return()
      }

      if(is.null(input$Process)){
        output$process_error <- renderText({
          validate(need(!is.null(input$Process), 'Please select at least 1 term (BP,MF,CC)'))
        })
        return()
      }
    
      isolate({
        df <- read_tsv(input$david_file$datapath)
        plot_info <- enrich_viz(df_d = df,
                                group = group,
                                CLUSTER_SCORE = input$CS,
                                FDR_SCORE = input$FDR,
                                GOTERM = input$Process,
                                CLUSTER = input$cluster)
        
        default_heatmap(plot_info[[1]],
                        row_font=input$font_r,
                        column_font=input$font_c,
                        clusterMethod=input$cluster,
                        g_mat = plot_info[[2]],
                        colorScale='Fixed',
                        legend_name='-log10(FDR)')
      })
    },height = input$height,width = input$width,)
  })
  
  # Saves parameters, tables, and figures in zip folder when user clicks download
  observeEvent(input$download, {
    v$doPlot <- input$plot
    
    output$error <- renderText({
      validate(need(v$doPlot == TRUE, 'Please generate initial plot before selecting download'))
    })
    if (v$doPlot == FALSE) return()
    
    ## create temporary directory
    zipdir <- tempfile("tmp", tmpdir = "./")
    if (!dir.exists(zipdir)) {
      dir.create(zipdir)  
    }
    setwd(zipdir)
    
    df <- read_tsv(input$david_file$datapath)
    plot_info <- enrich_viz(df_d = df,
                            group = group,
                            CLUSTER_SCORE = input$CS,
                            FDR_SCORE = input$FDR,
                            GOTERM = input$Process,
                            CLUSTER = input$cluster)
    
    # plot_info <- enrich_viz(df_d = df_d,
    #                         group = group,
    #                         CLUSTER = 'cj')

    label <- paste0(input$CS,"_FDRThreshold_",input$FDR,"_GOTERMS_",paste(input$Process,collapse = "_"),"_ClusterMethod_",input$cluster)

    ## Heatmap of DAVID enrichment with tsv
    if(!is.null(plot_info[[1]])){
      sig_mat <- as.data.frame(plot_info[[1]]) %>%
        mutate(GoTerms = rownames(x = .)) %>%
        relocate('GoTerms')
      write_tsv(x = sig_mat,file = paste0('FilteredSigEnrichments_NegLog10_ClusterScoreMin_',label,".tsv"))
      zip("data.zip",paste0('FilteredSigEnrichments_NegLog10_ClusterScoreMin_',label,".tsv"))

      ## Generate and save default heatmap
      pdf(file=paste0('HeatMap_ClusterScoreMin_',label,".pdf"),
          width = 16,height = 25)
      p1 <- default_heatmap(plot_info[[1]],
                            row_font=input$font_r,
                            column_font=input$font_c,
                            clusterMethod=input$cluster,
                            g_mat = plot_info[[2]],
                            colorScale='Fixed',
                            legend_name='-log10(FDR)')
      print(p1)
      dev.off()
      zip("data.zip",paste0('HeatMap_ClusterScoreMin_',label,".pdf"))
    }

    ## Gene Matrix tsv and heatmap
    if(!is.null(plot_info[[2]])){
      g_mat <- as.data.frame(plot_info[[2]]) %>%
        mutate(GoTerms = rownames(x = .)) %>%
        relocate('GoTerms')
      write_tsv(x = g_mat,file = paste0('GeneMatrixTable_ClusterScoreMin_',label,".tsv"))
      zip("data.zip",paste0('GeneMatrixTable_ClusterScoreMin_',label,".tsv"))

      pdf(file=paste0('GeneDistanceMap_ClusterScoreMin_',label,".pdf"),
          width = 30,height = 30)
      p4 <- gene_set_heatmap(1-plot_info[[2]])
      print(p4)
      dev.off()
      zip("data.zip",paste0('GeneDistanceMap_ClusterScoreMin_',label,".pdf"))
    }
    
    resulted_zip(paste0(zipdir,"/data.zip"))
    setwd(main_dir)
    
    runjs("$('#downloadData')[0].click();")
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("David_VisualizationAndSummary_", format(Sys.time(), "%a_%b_%d_%Y_%I-%M%p"), ".zip", sep="")
    },
    content = function(file) {
      file.copy(resulted_zip(), file)
    },
    contentType = "application/zip")
}

# Run the application 
shinyApp(ui = ui, server = server)
