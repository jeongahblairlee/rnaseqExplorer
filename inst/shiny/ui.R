ui <- fluidPage(
  theme = bs_theme(bootswatch = "lux"),


  titlePanel(tags$h2("Transcriptome Analysis Visualization", style = "color: #000000;
                     font-weight: bold;
                     font-size: 22px;")),

  tags$head(
    tags$style(HTML("
      .large-button {
        font-size: 20px !important;   /* Increase font size */
        padding: 10px 20px !important; /* Increase padding */
        background-color: #333333;     /* Charcoal black background color */
        color: white !important;        /* Button text color */
        border: none;                   /* Remove border */
        border-radius: 15px;            /* Rounded corners */
        cursor: pointer;                /* Pointer cursor on hover */
        box-shadow: 0 4px 8px rgba(0, 0, 0, 0.5); /* Subtle shadow */
        transition: background-color 0.3s, transform 0.3s; /* Smooth transition */
      }

      .large-button:hover {
        background-color: #555555;     /* Lighter shade on hover */
        transform: scale(1.05);        /* Slightly enlarge button on hover */
      }

      .large-button:active {
        transform: scale(0.95);        /* Slightly shrink button when clicked */
      }
    "))
  ),
  # Tabset for Navigation
  tabsetPanel(
    id = "tabs",
    type = "tabs",

    # Landing Page Tab
    tabPanel("Introduction",
             div(
               #style = "text-align: center; padding: 50px;",
               br(),
               uiOutput("intro_text"),  # Rendered output for the introduction text
               br(),  br(), br(), # Adds vertical space

               actionButton("go_to_analysis", "Go to Analysis",
                            class = "large-button")
             )
    ),

    tags$head(
      tags$style(HTML("
      .bold-text {
        font-weight: 800 !important;  /* Force bold text */
        color: black !important;       /* Change text color to black */
      }
          .large-text {
        font-size: 24px !important;    /* Increase font size */
      }
    "))
    ),

    # Analysis Page Tab
    tabPanel("Run analysis",
             fluidRow(
               column(
                 width = 4,
                 wellPanel(
                   style = "background-color: #ECF0F1;
                     border-radius: 30px;
                     padding: 15px;",
                   # Data Source Radio Buttons
                   radioGroupButtons("data_source", "Choose Data Source:",
                                     choices = list("Upload Data" = "upload", "Try Sample Data" = "sample"),
                                     selected = "upload",
                                     justified = TRUE,
                                     checkIcon = list(yes = icon("ok", lib = "glyphicon")),
                                     width = '100%',
                                     size = "sm"),  # Ensure it takes full width

                   tags$head(
                     tags$style(HTML("
      .custom-file-input .btn {
        padding: 5px 5px;      /* Adjust padding for button */
        width: 90px;
        height: 50px;
        font-size: 8pt;
        text-align: center;
        line-height: 50px;
      }
      .custom-file-input input[type='file'] {
        height: 50px;
        font-size: 5pt;
        text-align: center;
        line-height: 20px;
        width: 100%;            /* Full width for file input */
      }
    "))
                   ),

                   # File input (visible when "Upload Data" is selected)
                   div(class = "custom-file-input",
                       conditionalPanel(
                         condition = "input.data_source == 'upload'",
                         fileInput("file_input", "Upload CSV or Excel File:",
                                   width = '100%' ), # Full width
                         # helpText("Note: The file should contain gene expression data in the format of a count matrix.")
                       )),


                   radioGroupButtons(inputId = "plot_type",
                                     label = "Plot Type:",
                                     direction = "vertical",
                                     choices = list(
                                       "Violin plot" = "boxplot",
                                       "Heatmap" = "heatmap",
                                       "PCA" = "pca"
                                     ),
                                     selected = "boxplot",
                                     justified = FALSE,  # Ensure the buttons are not justified horizontally
                                     checkIcon = list(
                                       yes = icon("square-check"),
                                       no = icon("square")
                                     ),
                                     width = '100%',      # Full width
                                     size = "sm"         # Size can be adjusted as needed
                   ),

                   # Conditional UI for Box Plot
                   conditionalPanel(
                     condition = "input.plot_type == 'boxplot'",
                     radioGroupButtons("input_type", "Input Type:",
                                       choices = list("Select Name" = "select", "Enter Name" = "text"),
                                       selected = "select",
                                       justified = TRUE,
                                       width = '100%',
                                       size = "sm"),  # Full width

                     conditionalPanel(
                       condition = "input.input_type == 'select'",
                       selectizeInput("name", "Choose a name:",
                                      choices = NULL,
                                      width = '100%',
                                      size = "sm"),  # Full width
                       textOutput("select_warning")
                     ),


                     tags$head(
                       tags$style(HTML("
      .small-text-input .form-control {
        font-size: 10pt;    /* Set the desired font size */
      }
    "))
                     ),

                     conditionalPanel(
                       condition = "input.input_type == 'text'",
                       div(class = "small-text-input",
                           textInput("name",
                                     "Enter Name:",
                                     value = "",
                                     width = '100%')  # Full width
                       )
                     ),



                     selectInput("pvalue_type", "Choose p-value type:",
                                 choices = c("Adjusted p-value (DESeq2)" = "p.adj.deseq2", "Raw p-value (DESeq2)" = "p.deseq2")),  # Full width


                   ),

                   tags$head(
                     tags$style(HTML("

      .selectize-input {
        height: 50px;
        width: 100%;
        font-size: 12pt;
        padding-top: 5px;
      }

    "))
                   ),
                   tags$style(HTML("
    .selectize-input, .selectize-dropdown {
      font-size: 12px; /* Adjust this value to change the font size */
    }
  ")),


                   # Conditional UI for Heatmap
                   conditionalPanel(
                     condition = "input.plot_type == 'heatmap'",
                     div(class = "small-text-input",
                         numericInput("num_genes", "Choose number of top genes:",
                                      value = 10, min = 1, width = '100%')  # Full width
                     ),

                     radioGroupButtons("data_type", "Data Type:",
                                       choices = list("Ensembl ID" = "ensemblid", "Gene name" = "genename"),
                                       selected = "ensemblid", justified = TRUE,
                                       width = '100%',
                                       size = "sm"),  # Full width

                   ),



                   # Updated Run Analysis Button
                   actionBttn("run_analysis",
                              "Run Analysis",
                              style = "unite",  # Use 'fill' style for a cleaner look
                              color = "primary",  # Use primary color for better visibility
                              size = "md",  # Change to medium size for a balanced appearance
                              width = '100%')  # Full width
                 )
               ),

               # Output Plots
               column(
                 width = 8,
                 conditionalPanel(
                   condition = "input.plot_type == 'boxplot'",
                   plotOutput("boxPlot",
                              width = "75%",  # Ensure it takes full width
                              height = "400px"),
                   tableOutput("dataTable"),          # Download Button for Box Plot
                   downloadButton("download_boxplot", "Download Volin Plot",
                                  width = '100%',
                                  height = "50px",
                                  size = "md")
                 ),

                 conditionalPanel(
                   condition = "input.plot_type == 'heatmap'",
                   uiOutput("heatmap_ui"),
                   # Download Button for Heatmap
                   downloadButton("download_heatmap",
                                  "Download Heatmap Plot",
                                  width = '100%',
                                  height = "50px",
                                  size = "md"),
                   # Download Button for Heatmap matrix
                   downloadButton("download_matrix",
                                  "Download Heatmap Matrix",
                                  width = '100%',
                                  height = "50px",
                                  size = "md")
                 ),

                 conditionalPanel(
                   condition = "input.plot_type == 'pca'",
                   plotOutput("pcaPlot",
                              width = "100%",
                              height = "600px"),

                   # Download Button for PCA Plot
                   downloadButton("download_pca",
                                  "Download PCA Plot",
                                  width = '100%',
                                  height = "50px",
                                  size = "md")
                 )
               ),
             )
    )
  )
)
