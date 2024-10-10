

#' Launch Shiny Application
#'
#' This function initializes and runs the Shiny application by combining the UI and server components.
#' It serves as a wrapper around the `shinyApp()` function.
#' @param ui The UI layout of the Shiny application
#' @param server The server logic of the Shiny application
#' @examples
#' launchApp()
#' @export
# Define the function to launch the Shiny app
launchApp <- function() {
#   # List of required packages
#required_packages <- c(
#     "shiny", "shinyWidgets", "ggplot2", "dplyr", "rstatix", 
#     "ComplexHeatmap", "circlize", "DESeq2", "httr", "bslib"
#   )
#   
#   # Check if packages are installed, install if necessary
#   for (pkg in required_packages) {
#     if (!requireNamespace(pkg, quietly = TRUE)) {
#       install.packages(pkg)
#     }
#     library(pkg, character.only = TRUE)
#   }
#   


options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Setting a CRAN mirror

  # Define the required packages
  required_packages <- c(
    "shiny", "shinyWidgets", "ggplot2", 
    "dplyr", "rstatix", "circlize","bslib","readr"
  )
  
  # Check if packages are installed, install if necessary
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)  # Install from CRAN
    }
    library(pkg, character.only = TRUE)  # Load the package
  }
  
  if (!requireNamespace("BiocManager", quietly=TRUE)){
    install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
  }
  library(ComplexHeatmap)

  if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
   BiocManager::install("DESeq2")
  }
  library(DESeq2)


  if (!requireNamespace("httr", quietly=TRUE)){
  devtools::install_github("r-lib/httr")
  }
  library(httr)

 #    # Launch the shiny app
    shiny::shinyApp(ui = ui(), server = server)
  }
 #  
  
#' getting gene annotation
#'
#' exporting gene annotatioin based on ensembl id
#' @param ensembl_id 
#' @return gene annotation details
#' @examples 
#' temp1 <- get_gene_details('ENSG00000139618')
#' @export
get_gene_details <- function(ensembl_id) {
  base_url <- "https://rest.ensembl.org"
  endpoint <- paste0("/lookup/id/", ensembl_id, "?expand=1")

  response <- GET(paste0(base_url, endpoint), add_headers("Content-Type" = "application/json"))

  if (status_code(response) == 200) {
    content <- content(response, as = "parsed")
    return(content)
  } else {
    return(NULL)
  }
}


#' @export
citation.rnaseqExplorer <- function() {
  cat("If you use the rnaseqExplorer package in your research, please consider citing it as follows:\n")
  cat("Jeongah Lee. (2024). rnaseqExplorer: A comprehensive Shiny application for analyzing RNA sequencing data. R package version 0.1.0. https://github.com/jeongahblairlee/rnaseqExplorer\n")
}

#' @export
citation <- function(package = "rnaseqExplorer") {
  if (package == "rnaseqExplorer") {
    citation.rnaseqExplorer()
  } else {
    methods::citation(package)
  }
}

.onLoad <- function(libname, pkgname) {
  print("Package rnaseqExplorer has been loaded!")  # Add this line for debugging
  citation.rnaseqExplorer()
}





#' Creating the Shiny UI for displaying user information
#'
#' Generates a Shiny UI layout for displaying user details such as name, email, and activities.
#' This function defines the structure and layout of the app's UI.
#' @examples
#' ui <- create_user_ui()
#' @export
ui <- function() {
  shiny::fluidPage(
    theme = bs_theme(bootswatch = "lux"),
    
    # Title Panel
    titlePanel("rnaseqExplorer"),
    
    # Add CSS for consistent spacing and styles
    tags$head(
      tags$style(HTML("
        .large-button {
          font-size: 20px !important;
          padding: 10px 20px !important;
          background-color: #333333;
          color: white !important;
          border: none;
          border-radius: 15px;
          cursor: pointer;
          box-shadow: 0 4px 8px rgba(0, 0, 0, 0.5);
          transition: background-color 0.3s, transform 0.3s;
        }
        
        .large-button:hover {
          background-color: #555555;
          transform: scale(1.05);
        }

        .large-button:active {
          transform: scale(0.95);
        }

        .tab-panel {
          margin-bottom: 20px;
        }

        .bold-text {
          font-weight: 800 !important;
          color: black !important;
        }

        .large-text {
          font-size: 24px !important;
        }
      "))
    ),
    
    # Tabset for Navigation
    tabsetPanel(
      id = "tabs",
      type = "tabs",
      
      # Landing Page Tab
      tabPanel("Introduction",
               div(class = "tab-panel",
                   br(),
                   uiOutput("intro_text"),  # Rendered output for the introduction text
                   br(), br(), br(),
                   actionButton("go_to_analysis", "Go to Analysis", class = "large-button")
               )
      ),
      
      # Workflow Tab
      tabPanel("Workflow",
               div(class = "tab-panel",
                   br(),
                   h4("Workflow"),
                   br(),
                   tags$img(src = "https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/workflow.png",
                            style = "width: 100%; height: auto; display: block; margin: auto;",
                            alt = "workflow"),
                   br()
               )
      ),
      
      # How to Use Tab
      tabPanel("How to use",
               div(class = "tab-panel",
                   br(),
                   h4("1. Upload your data:"),
                   tags$img(src = "https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/function1.png",
                            style = "width: 100%; height: auto; display: block; margin: auto;",
                            alt = "function1"),
                   br(),
                   br(),
                   br(),
                   h4("2. Explore the data:"),
                   h5("2-1. Violin plot:"),
                   tags$img(src = "https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/function2.png",
                            style = "width: 100%; height: auto; display: block; margin: auto;",
                            alt = "function2"),
                   br(),
                   br(),
                   br(),
                   h5("2-2. Heatmap:"),
                   tags$img(src = "https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/function3.png",
                            style = "width: 100%; height: auto; display: block; margin: auto;",
                            alt = "function3"),
                   br(),
                   br(),
                   br(),
                   h5("2-3. PCA:"),
                   tags$img(src = "https://raw.githubusercontent.com/jeongahblairlee/rnaseqExplorer/refs/heads/main/notebook/function4.png",
                            style = "width: 100%; height: auto; display: block; margin: auto;",
                            alt = "function4"),
                   br()
               )
      ),
      
      # Analysis Page Tab
      tabPanel("Run analysis",
               div(class = "tab-panel",
                   fluidRow(
                     column(
                       width = 4,
                       wellPanel(
                         style = "background-color: #ECF0F1; border-radius: 30px; padding: 15px;",
                         # Data Source Radio Buttons
                         radioGroupButtons("data_source", "Choose Data Source:",
                                           choices = list("Upload Data" = "upload", "Try Sample Data" = "sample"),
                                           selected = "upload",
                                           justified = TRUE,
                                           checkIcon = list(yes = icon("ok", lib = "glyphicon")),
                                           width = '100%',
                                           size = "sm"),
                         
                         # File input (visible when "Upload Data" is selected)
                         div(class = "custom-file-input",
                             conditionalPanel(
                               condition = "input.data_source == 'upload'",
                               fileInput("file_input", "Upload CSV or Excel File:", width = '100%') # Full width
                             )),
                         
                         radioGroupButtons(inputId = "plot_type", label = "Plot Type:",
                                           direction = "vertical",
                                           choices = list(
                                             "Violin plot" = "boxplot",
                                             "Heatmap" = "heatmap",
                                             "PCA" = "pca"
                                           ),
                                           selected = "boxplot",
                                           justified = FALSE,
                                           checkIcon = list(yes = icon("square-check"), no = icon("square")),
                                           width = '100%',
                                           size = "sm"),
                         
                         # Conditional UI for Box Plot
                         conditionalPanel(
                           condition = "input.plot_type == 'boxplot'",
                           radioGroupButtons("input_type", "Input Type:",
                                             choices = list("Select Name" = "select", "Enter Name" = "text"),
                                             selected = "select",
                                             justified = TRUE,
                                             width = '100%',
                                             size = "sm"),
                           
                           conditionalPanel(
                             condition = "input.input_type == 'select'",
                             selectizeInput("name", "Choose a name:", choices = NULL, width = '100%', size = "sm"),
                             textOutput("select_warning")
                           ),
                           
                           conditionalPanel(
                             condition = "input.input_type == 'text'",
                             div(class = "small-text-input",
                                 textInput("name", "Enter Name:", value = "", width = '100%'))
                           ),
                           
                           selectInput("pvalue_type", "Choose p-value type:",
                                       choices = c("Adjusted p-value (DESeq2)" = "p.adj.deseq2", "Raw p-value (DESeq2)" = "p.deseq2"))
                         ),
                         
                         # Conditional UI for Heatmap
                         conditionalPanel(
                           condition = "input.plot_type == 'heatmap'",
                           div(class = "small-text-input",
                               numericInput("num_genes", "Choose number of top genes:", value = 10, min = 1, width = '100%')),
                           
                           radioGroupButtons("data_type", "Data Type:",
                                             choices = list("Ensembl ID" = "ensemblid", "Gene name" = "genename"),
                                             selected = "ensemblid", justified = TRUE,
                                             width = '100%',
                                             size = "sm")
                         ),
                         
                         # Updated Run Analysis Button
                         actionBttn("run_analysis", "Run Analysis",
                                    style = "unite",
                                    color = "primary",
                                    size = "md",
                                    width = '100%')
                       )
                     ),
                     
                     # Output Plots
                     column(
                       width = 8,
                       conditionalPanel(
                         condition = "input.plot_type == 'boxplot'",
                         plotOutput("boxPlot", width = "450px", height = "450px"),
                         tableOutput("dataTable"),          # Download Button for Box Plot
                         downloadButton("download_boxplot", "Download Violin Plot", width = '100%', height = "50px", size = "md")
                       ),
                       
                       conditionalPanel(
                         condition = "input.plot_type == 'heatmap'",
                         uiOutput("heatmap_ui"),
                         # Download Button for Heatmap
                         downloadButton("download_heatmap", "Download Heatmap Plot", width = '100%', height = "50px", size = "md"),
                         # Download Button for Heatmap matrix
                         downloadButton("download_matrix", "Download Heatmap Matrix", width = '100%', height = "50px", size = "md")
                       ),
                       
                       conditionalPanel(
                         condition = "input.plot_type == 'pca'",
                         plotOutput("pcaPlot", width = "700px", height = "600px"),
                         
                         # Download Button for PCA Plot
                         downloadButton("download_pca", "Download PCA Plot", width = '100%', height = "50px", size = "md")
                       )
                     )
                   )
               )
      )
    )
  )
}








# You can create your server logic here if needed



#' Defining server logic for handling user interactions and data
#'
#' Implements the server-side logic, including reactive expressions and observers,
#' to process and display user information in the Shiny app.
#' @examples
#' server <- function(input, output, session) {
#'   server_logic()
#' }
#' @export
server <- function(input, output, session) {

  values <- reactiveValues(processed_data = NULL)

  # Set max upload size to 100MB
  options(shiny.maxRequestSize = 100*1024^2)

  # Reactive function to load either sample (airway) or uploaded data
  dataset <- reactive({
    req(input$data_source)  # Ensure input is available

    if (input$data_source == "sample") {
      # Use the airway dataset (sample data)
      withProgress(message = "Processing data...", value = 0, {

        # Stage 1: Loading airway dataset
        incProgress(0.1, detail = "Loading sample dataset...")
        library("airway")
        data(airway)

        # Stage 2: Processing count data
        incProgress(0.2, detail = "Processing count data...")
        cts <- assay(airway)
        cts <- cts[1:500,]
        col_names <- colnames(cts)
        ordered_col_names <- col_names[order(col_names)]
        cts <- cts[, ordered_col_names]

        # Stage 3: Processing metadata
        incProgress(0.1, detail = "Processing metadata...")
        coldata <- colData(airway)
        coldata <- coldata[order(rownames(coldata)), ]
        colnames(coldata)[3] <- "condition"
        coldata$rowname <- rownames(coldata)
        coldata$condition <- ifelse(coldata$condition == "untrt", "control", "treatment")
        coldata$condition <- factor(coldata$condition, levels = c("control", "treatment"))

        list(cts = cts, coldata = coldata)
      })

    } else if (!is.null(input$file_input)) {
      # Use user-uploaded data (CSV or Excel)
      withProgress(message = "Uploading and processing file...", value = 0, {

        file <- input$file_input
        ext <- tools::file_ext(file$datapath)

        tryCatch({
          # Stage 1: Reading file
          incProgress(0.1, detail = "Reading file...")
          if (ext == "csv") {
            data <- read.csv(file$datapath)
          } else if (ext == "xlsx") {
            data <- readxl::read_excel(file$datapath, 1)  # Read the first sheet
          } else {
            stop("Unsupported file format.")  # Throw an error for unsupported formats
          }

          # Stage 2: Converting data to data frame and removing duplicates
          incProgress(0.2, detail = "Processing data frame...")
          df <- as.data.frame(data)
          df <- df[!duplicated(df[, 1]), ]

          # Stage 3: Setting up row names and metadata
          incProgress(0.2, detail = "Formatting metadata...")
          rownames(df) <- df[, 1]
          rownames(df)[1] <- "condition"
          df[, 1] <- NULL
          coldata <- t(df[1, ])

          # Stage 4: Processing count matrix
          incProgress(0.2, detail = "Processing count matrix...")
          cts <- as.data.frame(df[-1, ])
          cts <- mutate_all(cts, function(x) as.integer(as.character(x)))
          cts <- na.omit(cts)

          # Remove rows with all 0 values
          incProgress(0.1, detail = "Removing rows with all zero values...")
          cts <- cts[rowSums(cts) != 0, ]

          # Stage 5: Finalizing column names and matrix validation
          incProgress(0.2, detail = "Validating data...")
          colnames(coldata) <- "condition"
          if (!all(rownames(coldata) %in% colnames(cts))) {
            stop("Error: Column names in the metadata do not match the count matrix.")  # Format error message
          }

          cts <- cts[, rownames(coldata)]
          list(cts = cts, coldata = coldata)

        }, error = function(e) {
          # Handle errors and display warning message
          showNotification(paste("Error loading dataset:", e$message,
                                 "Please check the dataset format."), type = "error")
          return(NULL)
        })
      })
    }
  })



  observeEvent(input$run_analysis, {
    data <- dataset()

    # Validate that dataset and the necessary variables exist
    validate(
      need(!is.null(data), "Error: Dataset is missing!"),
      need("cts" %in% names(data), "Error: Count data (cts) is missing in the dataset!"),
      need("coldata" %in% names(data), "Error: Column data (coldata) is missing in the dataset!")
    )

    cts <- data$cts
    coldata <- data$coldata

    withProgress(message = "Processing data...", value = 0, {
      incProgress(0.3, detail = "Running DESeq2 analysis...")

      # Ensure cts and coldata are not NULL or empty
      validate(
        need(!is.null(cts) && nrow(cts) > 0, "Error: cts is empty!"),
        need(!is.null(coldata) && nrow(coldata) > 0, "Error: coldata is empty!")
      )

      # Create DESeq dataset and run DESeq2 analysis
      dds <- tryCatch({
        DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)
      }, error = function(e) {
        showNotification(paste("Error in DESeqDataSet creation:", e$message), type = "error")
        return(NULL)
      })


      # If dds is NULL due to error, stop processing
      validate(need(!is.null(dds), "Error: Failed to create DESeq dataset!"))

      dds <- DESeq(dds)

      incProgress(0.2, detail = "Filtering and normalizing data...")

      min_samples <- ceiling(ncol(dds) / 2)
      keep <- rowSums(counts(dds) >= 10) >= min_samples
      dds_filtered <- dds[keep, ]


      # Check if filtering resulted in valid dds_filtered
      validate(
        need(nrow(dds_filtered) > 0, "Error: No rows left after filtering!")
      )

      res <- results(dds_filtered)
      resOrdered <- res[order(res$padj),]


      # Run vst transformation alone to see if it throws an error
      vsd <- tryCatch({
        # Attempt variance stabilizing transformation (VST)
        assay(vst(dds_filtered))
      }, error = function(e) {
        # Print the error message
        message(paste("Switching to rlog due to small sample size.", e$message))

        # Fall back to rlog
        rlog_assay <- tryCatch({
          assay(rlog(dds_filtered))
        }, error = function(r_error) {
          message(paste("Error in rlog:", r_error$message))
          return(NULL)  # Return NULL if rlog also fails
        })

        return(rlog_assay)  # Return the rlog-transformed assay
      })



      # If vsd is NULL due to error, stop processing
      validate(need(!is.null(vsd), "Error: Failed to perform variance stabilizing transformation!"))

      Z <- t(scale(t(vsd)))

      incProgress(0.2, detail = "Reshaping data for plotting...")

      # Reshape data for plotting and check for errors
      library(reshape)
      m <- tryCatch({
        melt(vsd)
      }, error = function(e) {
        showNotification("Error in reshaping data for plotting: " , e$message, type = "error")
        return(NULL)
      })

      validate(need(!is.null(m), "Error: Failed to reshape data!"))

      colnames(m) <- c("gene", "rowname", "value")
      coldata <- as.data.frame(coldata)
      coldata$rowname <- rownames(coldata)

      m2 <- merge(m, coldata, by = "rowname")
      incProgress(0.1, detail = "Wrapping up data processing...")

      # Return the processed data as a list
      values$processed_data<- list(
        coldata = coldata,
        cts = cts,
        resOrdered = resOrdered,
        res = res,  # Ensure res is returned for access in other plots
        dds_filtered = dds_filtered,
        vsd = vsd,
        m = m,
        Z = Z,
        m2 = m2
      )
    })
  })




  # Reset processed_data when input data changes
  observeEvent({
    input$data_source
    input$file_input
  }, {
    values$processed_data <- NULL
  })

  observe({
    if (!is.null(values$processed_data)) {
      print("Processed data has been updated.")
      print(head(values$processed_data$resOrdered))  # Check the contents of resOrdered
    }
  })

  # Update selectizeInput choices dynamically
  observe({
    data <- values$processed_data

    validate(
      need(!is.null(data$resOrdered), "No processed data available yet.")
    )
    top_n <- 10000000  # Limit the number of choices to display
    top_genes <- rownames(data$resOrdered)

    if (length(top_genes) > top_n) {
      top_genes <- head(top_genes, top_n)
      output$select_warning <- renderText({
        "Warning: The select input contains a large number of options;
        please consider using the search functionality to find your gene of interest."
      })
    } else {
      output$select_warning <- renderText({ "" })
    }
    # Update selectizeInput with the list of top_genes
    updateSelectizeInput(session, "name",
                         choices = top_genes,
                         selected = NULL,  # Default selection is NULL
                         server = TRUE)

  })

  output$boxPlot <- renderPlot({
    req(input$plot_type == "boxplot")

    data <- values$processed_data

    validate(
      need(!is.null(data$m2), "No processed data available yet.")
    )
    selected_name <- input$name
    m2 <- data$m2
    res <- data$res
    # Filter the dataset based on selected name
    m3 <- m2[m2$gene == selected_name, ]
    # Check if the dataset is not empty
    if (nrow(m3) < 2) {
      return(NULL)
    }
    # m3$condition <- factor(m3$condition, levels <- c("untreated", "treated"))
    library(ggpubr)
    pd <- position_dodge(width = 0.5)
    p <- ggviolin(m3, x = "condition", y = "value", fill = "condition", alpha = 0.75, size = 0.7) +
      ylab("Gene expression level") +
      stat_boxplot(geom = 'errorbar', position = pd, linetype = 1, width = 0.2) +
      geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.7, position = pd, fill = "black")

    library(rstatix)
    m3 <- as.data.frame(m3)
    stat.test <- m3 %>%
      t_test(value ~ condition) %>%
      adjust_pvalue() %>%
      add_significance("p.adj") %>%
      add_xy_position(x = "condition")

    stat.test$p.adj.deseq2 <- round(res[rownames(res) == selected_name, 6], digits = 5)
    stat.test$p.deseq2 <- round(res[rownames(res) == selected_name, 5], digits = 5)
    stat.test$p.adj.deseq2_full <- res[rownames(res) == selected_name, 6]
    stat.test$p.deseq2_full <- res[rownames(res) == selected_name, 5]

    pvalue_label <- input$pvalue_type
    warning_message <- ""

    if (pvalue_label == "p.adj.deseq2") {
      if (is.na(stat.test$p.adj.deseq2)) {
        warning_message <- "The calculation of adjusted p-value is not possible."
      }
      pp <- p + stat_pvalue_manual(stat.test,
                                   label = pvalue_label,
                                   size = 4,
                                   bracket.shorten = 0.1) +
        scale_fill_manual(values = c("gold1", "darkmagenta")) +
        labs(title = paste(m3$gene[1]),
             caption = bquote(italic(.(paste0("adjusted P-value: ",stat.test$p.adj.deseq2_full)))),
             fill = "Condition") +
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 12),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = 14, face = "bold"),
              axis.title.y = element_text(size = 15),
              legend.text = element_text(size = 12),
              legend.position = "right",
              legend.title = element_text(size = 12, face = "bold"),
              plot.title = element_text(size = 15, face = "bold",hjust = 0.5, vjust = 2),
              plot.caption = element_text(size = 12, vjust = 0.1, family="italic"))
    } else if (pvalue_label == "p.deseq2") {
      if (is.na(stat.test$p.deseq2)) {
        warning_message <- "The calculation of raw p-value is not possible."
      }
      pp <- p + stat_pvalue_manual(stat.test,
                                   label = pvalue_label,
                                   size = 4,
                                   bracket.shorten = 0.1) +
        scale_fill_manual(values = c("gold1", "darkmagenta")) +
        labs(title = paste(m3$gene[1]),
             caption = bquote(italic(.(paste0("P-value: ",stat.test$p.deseq2_full)))),
             fill = "Condition") +
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 12),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size = 14, face = "bold"),
              axis.title.y = element_text(size = 15),
              legend.text = element_text(size = 12),
              legend.position = "right",
              legend.title = element_text(size = 12, face = "bold"),
              plot.title = element_text(size = 15, face = "bold",hjust = 0.5,vjust = 2),
              plot.caption = element_text( size = 12, vjust = 0.1, family= "italic"))
    }

    # if (warning_message != "") {
    #   output$warning <- renderText({
    #     warning_message
    #   })
    print(pp)
    # } else {
    #   output$warning <- renderText({
    #     ""
    #   })
    #   print(pp)
    # }

    # Download Boxplot
    output$download_boxplot <- downloadHandler(
      filename = function() {
        paste("boxplot", "_",m3$gene[1],"_",Sys.Date(), ".tiff", sep = "")
      },
      content = function(file) {
        tiff(file,res=300,units="in", width=5.5, height=5.5)
        print(pp)
        dev.off()
      }
    )

  })


  # Render the gene details as a table
  output$dataTable <- renderTable({
    req(input$plot_type == "boxplot")
    data <- values$processed_data
    selected_name <- input$name

    # Fetch gene details
    gene_details <- get_gene_details(selected_name)

    # Check if gene_details is not NULL and contains the necessary fields
    if (!is.null(gene_details)) {
      # Safely access each field with ifelse to avoid errors if a field is missing
      gene_info <- data.frame(
        "Ensembl_ID"  = ifelse(!is.null(gene_details$id), gene_details$id, NA),
        "Gene_name"   = ifelse(!is.null(gene_details$display_name), gene_details$display_name, NA),
        "Chr"         = ifelse(!is.null(gene_details$seq_region_name), gene_details$seq_region_name, NA),
        #"Start"      = ifelse(!is.null(gene_details$start), gene_details$start, NA),  # Can add back if needed
        #"End"        = ifelse(!is.null(gene_details$end), gene_details$end, NA),      # Can add back if needed
        "Description" = ifelse(!is.null(gene_details$description), gene_details$description, NA),
        #"Species"    = ifelse(!is.null(gene_details$species), gene_details$species, NA),  # Can add back if needed
        #"Assembly"   = ifelse(!is.null(gene_details$assembly_name), gene_details$assembly_name, NA),  # Can add back if needed
        stringsAsFactors = FALSE
      )

    } else {
      gene_info <- data.frame("Error" = "Gene details not found")
    }

    gene_info
  })




  # Dynamically set heatmap plot dimensions based on top genes selected
  output$heatmap_ui <- renderUI({
    req(input$plot_type == "heatmap")

    # Calculate height dynamically based on the number of genes selected
    num_genes <- as.numeric(input$num_genes)
    data <- values$processed_data
    Z <- data$Z
    num_cols <- ncol(Z)
    plot_width <- paste0(500 + num_cols * 10, "px")
    plot_height <- paste0(300 + num_genes * 10, "px")  # Adjust height scaling factor as necessary

    # Return the plot output with dynamic height
    plotOutput("heatmapPlot", width = plot_width, height = plot_height)
  })



  # New Heatmap Plot code
  output$heatmapPlot <- renderPlot({
    req(input$plot_type == "heatmap")
    data <- values$processed_data
    res <- data$res
    Z <- data$Z
    coldata <- data$coldata
    # Getting the number of top genes based on user input
    num_genes <- as.numeric(input$num_genes)

    withProgress(message = 'Processing data for heatmap...', {

      # Filter out rows with NA p-values
      res_filtered <- res[!is.na(res$pvalue), ]

      # Sort by p-value in ascending order
      res_sorted <- res_filtered[order(res_filtered$pvalue), ]

      # Select the top genes
      top10_genes <- head(res_sorted, num_genes)
      top10_cts <- Z[rownames(Z) %in% rownames(top10_genes), ]

      col_fun = colorRamp2(seq(min(top10_cts), max(top10_cts), length = 3),
                           c("royalblue3", "#EEEEEE", "firebrick2"), space = "RGB")

      # Extract unique condition levels from coldata
      unique_conditions <- unique(coldata$condition)

      # Define colors for each unique condition
      condition_colors <- setNames(c("gold1", "darkmagenta"), unique_conditions)

      # Ensure that coldata$condition has the same levels as the names of condition_colors
      anno_col <- list(Condition = condition_colors)

      ha = HeatmapAnnotation(Condition = coldata$condition,
                             col = anno_col,
                             show_annotation_name = TRUE)

      data_type <- input$data_type
      if (data_type == "genename") {
        names <- rownames(top10_cts)
        gene_symbol <-c()
        withProgress(message = 'Fetching gene names...', {
          for (i in 1:length(names)) {
            # incProgress(1 / length(names), detail = paste("Fetching:", names[i]))  # Update progress
            gene_name <- get_gene_details(names[i])

            if (!is.null(gene_name) && !is.null(gene_name$id) && !is.null(gene_name$display_name)) {
              # If both id and display_name are not NULL, assign them
              gene_symbol$id[i] <- gene_name$id
              gene_symbol$symbol[i] <- gene_name$display_name
            } else {
              # If gene_name is NULL or id/display_name are NULL, assign id to both
              gene_symbol$id[i] <- names[i]  # Use the current gene name as the id
              gene_symbol$symbol[i] <- gene_symbol$id[i]  # Set symbol to id
              }
            }
        gene_symbol <- as.data.frame(gene_symbol)
        })


        
        # Extract relevant columns from gene_symbol
        colnames(gene_symbol) <- c("id", "display_name")

        # Convert to a named vector for easier lookup
        display_name_map <- setNames(gene_symbol$display_name, gene_symbol$id)

        # Replace row names in top10_cts
        rownames(top10_cts) <- display_name_map[rownames(top10_cts)]
        heatmap_plot <- Heatmap(top10_cts, name = "Z-score",
                                column_title = paste0("TOP ", num_genes, " gene expression"),
                                col = col_fun,
                                row_names_gp = gpar(fontsize = 10),
                                show_row_names = TRUE,
                                top_annotation = ha,
                                show_column_names = TRUE,
                                show_row_dend = FALSE,
                                show_column_dend = TRUE,
                                column_names_gp = gpar(fontsize = 10),
                                row_names_side = "left")

        draw(heatmap_plot, merge_legend = TRUE, heatmap_legend_side = "right",
             annotation_legend_side = "right")
      }else{

        heatmap_plot <- Heatmap(top10_cts, name = "Z-score",
                                column_title = paste0("TOP ", num_genes, " gene expression"),
                                col = col_fun,
                                row_names_gp = gpar(fontsize = 10),
                                show_row_names = TRUE,
                                top_annotation = ha,
                                show_column_names = TRUE,
                                show_row_dend = FALSE,
                                show_column_dend = TRUE,
                                column_names_gp = gpar(fontsize = 10),
                                row_names_side = "left")

        draw(heatmap_plot, merge_legend = TRUE, heatmap_legend_side = "right",
             annotation_legend_side = "right")
      }
    })

    # Download Heatmap
    output$download_heatmap <- downloadHandler(
      filename = function() {
        paste("heatmap_top", num_genes,"_genes_", Sys.Date(), ".tiff", sep = "")
      },
      content = function(file) {
        num_genes <- as.numeric(input$num_genes)
        data <- values$processed_data
        Z <- data$Z
        num_cols <- ncol(Z)
        plot_w <- paste0((500 + num_cols * 10) / 96)
        plot_h <- paste0((300 + num_genes * 10) / 96)
        tiff(file, res = 300, units = "in", width = plot_w, height = plot_h)
        heatmap_plot <- Heatmap(top10_cts, name = "Z-score",
                                column_title = paste0("TOP ", num_genes, " gene expression"),
                                col = col_fun,
                                row_names_gp = gpar(fontsize = 10),
                                show_row_names = TRUE,
                                top_annotation = ha,
                                show_column_names = TRUE,
                                show_row_dend = FALSE,
                                show_column_dend = TRUE,
                                column_names_gp = gpar(fontsize = 10),
                                row_names_side = "left")
        draw(heatmap_plot, merge_legend = TRUE, heatmap_legend_side = "right",
             annotation_legend_side = "right")
        dev.off()
      }
    )

  })

  # Download Heatmap matrix
  output$download_matrix <- downloadHandler(

    filename = function() {
      paste("heatmap_top", input$num_genes, "_gene_expression_", Sys.Date(), ".csv", sep = "")
    },

    content = function(file) {
      # Retrieve number of genes
      num_genes <- as.numeric(input$num_genes)
      # Ensure 'res' and 'vsd' are defined in the global environment or accessible
      data <- values$processed_data
      res <- data$res
      vsd <- data$vsd
      #bring the res
      res_filtered <- res[!is.na(res$pvalue), ]
      
      # Sort by p-value in ascending order
      res_sorted <- res_filtered[order(res_filtered$pvalue), ]
      
      # Select the top genes
      top10_genes <- head(res_sorted, num_genes)
      top10_vsd <- vsd[rownames(vsd) %in% rownames(top10_genes), ]
      
      names <- rownames(top10_vsd)
      
      # Initialize an empty data frame with the appropriate column names
      gene_symbol <- data.frame(id = character(), display_name = character(), stringsAsFactors = FALSE)
      
      # Loop over the 'names' vector
      for (i in seq_along(names)) {
        
        # Retrieve gene details for the ith gene
        gene_name <- get_gene_details(names[i])
        
        # If gene details are not null, append the values to gene_symbol
        if (!is.null(gene_name)) {
          gene_symbol <- rbind(gene_symbol,
                               data.frame(id = ifelse(!is.null(gene_name$id), gene_name$id, NA),
                                          display_name = ifelse(!is.null(gene_name$display_name), gene_name$display_name, NA),
                                          stringsAsFactors = FALSE))
        }
      }
      
      
      # Extract relevant columns from gene_symbol
      colnames(gene_symbol) <- c("id", "gene_name")
      top10_vsd <- as.data.frame(top10_vsd)
      
      # Create a named vector for display_name
      display_name_vector <- setNames(gene_symbol$gene_name, gene_symbol$id)
      
      # Add the display_name to top10_vsd
      top10_vsd$gene_name <- display_name_vector[row.names(top10_vsd)]
      
      # Convert the DESeqResults object to a data frame
      res_df <- as.data.frame(res_sorted)
      res_df <- res_df[, c("log2FoldChange", "pvalue", "padj")]
      
      # Add row names as a column to the res_df for merging
      res_df$ID <- row.names(res_df)
      
      # Add row names as a column to top10_vsd for merging
      top10_vsd$ID <- row.names(top10_vsd)
      
      # Perform the merge using left_join
      final_result <- top10_vsd %>%
        dplyr::left_join(res_df, by = "ID")
      
      # Clean up the final result
      rownames(final_result) <- final_result$ID
      final_result$ID <- NULL
      
      # Rearrange columns
      final_result <- final_result %>%
        dplyr::select(gene_name, everything())
      
      # Write the final result to a CSV file
      write_csv(final_result, file)
    }
  )
  


  output$pcaPlot <- renderPlot({
    req(input$plot_type == "pca")

    # Initialize the progress bar
    withProgress(message = 'Generating PCA plot...', value = 0, {

      # Increment progress (e.g., start at 0.1)
      incProgress(0.1, detail = "Loading dataset...")
      data <- values$processed_data
      vsd <- data$vsd
      coldata <- data$coldata
      # Increment progress after loading the data
      incProgress(0.2, detail = "Performing PCA analysis...")
      library(PCAtools)

      p <- pca(vsd, scale = FALSE, metadata = coldata, removeVar = 0.1)

      unique_conditions <- unique(coldata$condition)

      # Define colors for each unique condition
      condition_colors <- setNames(
        c("gold1", "darkmagenta"),  # List colors corresponding to the conditions
        unique_conditions          # Names of the colors
      )

      # Increment progress after PCA calculation
      incProgress(0.9, detail = "Generating PCA plot...")

      pca_plot <- biplot(p, lab = p$yvars, labSize = 4,
                         title = 'Principal component analysis',
                         colby = 'condition',
                         colLegendTitle = "Condition",
                         legendPosition = 'right', legendLabSize = 12, legendIconSize = 7,
                         colkey = condition_colors,
                         sizeLoadingsNames = 3.5,
                         colLoadingsNames = "darkblue",
                         boxedLoadingsNames = TRUE,
                         showLoadings = FALSE,
                         showLoadingsNames = TRUE,
                         pointSize = 5.5,
                         max.overlaps = 10,
                         ellipse = FALSE,
                         ellipseLevel = 0.9,
                         ellipseFill = TRUE,
                         ellipseAlpha = 0.2,
                         ellipseLineSize = 0.00001)

      # Final stage of progress
      incProgress(0.8, detail = "Finalizing plot...")

      # Print PCA plot
      print(pca_plot)

      # Finish progress
      incProgress(1, detail = "Done")
    })


    # Download PCA Plot
    output$download_pca <- downloadHandler(
      filename = function() {
        paste("pca_plot", "_",Sys.Date(), ".tiff", sep = "")
      },
      content = function(file) {
        tiff(file,res=300,units="in", width=7.5, height=7)

        pca_plot <- biplot(p, lab = p$yvars, labSize = 4,
                           title = 'Principal component analysis',
                           colby = 'condition',
                           colLegendTitle = "Condition",
                           legendPosition = 'right', legendLabSize = 14, legendIconSize = 10,
                           colkey = condition_colors,
                           sizeLoadingsNames = 3.5,
                           colLoadingsNames = "darkblue",
                           boxedLoadingsNames = TRUE,
                           showLoadings = FALSE,
                           showLoadingsNames = TRUE,
                           pointSize = 5.5,
                           max.overlaps = 10,
                           ellipse = FALSE,
                           ellipseLevel = 0.9,
                           ellipseFill = TRUE,
                           ellipseAlpha = 0.2,
                           ellipseLineSize = 0.00001)


        # Print PCA plot
        print(pca_plot)
        dev.off()
      }
    )

  })



  output$intro_text <- renderText({
    HTML(paste0(
      
      "<span class='large-text'>Welcome to <span class='bold-text'>rnaseqExplorer</span> </span>",
      br(), br(),
      "<span class='bold-text'>rnaseqExplorer</span> is a Bioconductor package that includes a Shiny application for analyzing gene expression data. 📊<br>",
      br(),
      "This tool serves as a general-purpose interactive visualization for RNA-seq analysis, guiding users in easily exploring transcriptome data. 🌐<br>",
      br(),
      "<span class='bold-text'>rnaseqExplorer</span> provides functionalities for data analysis, including filtering, normalization, and understanding gene expression through visualizations like:<br>",
      "<ul>",  # Start an unordered list
      "<li>🎻 Violin plots</li>",
      "<li>🗺️ Heatmaps</li>",
      "<li>📈 PCA</li>",
      "<li>🔍 Gene annotation</li>",
      "</ul>",  # End the unordered list
      br(),
      "With its dynamic interactive design, <span class='bold-text'>rnaseqExplorer</span> is an essential tool for RNA-seq dataset analysis, empowering bench biologists to conduct exploratory data analysis with ease while delivering in-depth insights for seasoned data analysts. 💡"
    ))
  })

  # Observe event for switching tabs
  observeEvent(input$go_to_analysis, {
    updateTabsetPanel(session, "tabs", selected = "Run analysis")
  })

}

shinyApp(ui = ui, server = server)



