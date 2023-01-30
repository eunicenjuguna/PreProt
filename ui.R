
#install.packages(c("shiny","tidyverse","DT","shinydashboard","shinydashboardPlus"))

library(shiny)
library(tidyverse)
library(DT)
library(shinydashboard)
library(shinydashboardPlus)

proteomics_contaminats<- c("REV","CON","amalyse","keratin","abumin","trypsin","lysin")

header <- dashboardHeader(
  title = "PreProt",
  titleWidth = 300
)

sidebar <- dashboardSidebar(
  width = 300,
  menuItem("TMT Data Preparation", tabName = "tmt_preparation",
           menuSubItem("Upload TMT Raw Data", 
                       tabName = "tmt_data_upload_tab"),
           menuSubItem("Filter Preferred Intensity Columns",
                       tabName = "tmt_reportercorrected_tab")
  ),
 
  menuItem("LFQ Data Preparation", tabName = "lfq_analysis",
           menuSubItem("Upload LFQ Raw Data",
                       tabName = "lfq_data_upload_tab"),
           menuSubItem("Filter LFQ Intensity Columns",
                       tabName = "lfq_intensity_tab")
           ),
  br(),
  hr(),
  
  radioButtons(inputId = "select_analysis_type", label = "choose analysis type",choices = c("TMT","LFQ")),
  
  br(),
  hr(),
  
  menuItem(" Data Cleaning", tabName = "data cleaning",
           menuSubItem("Filter by Threshold",
                       tabName = "filter_missing_tab"),
           menuSubItem("Filter the Contaminants ",
                       tabName = "contaminants_filter_tab")
  ),
  menuItem("Imputation", tabName = "imputation",
           menuSubItem("Inspect Missing Values",
                       tabName = "mv_inspection_tab"),
           menuSubItem("Check Missing Value Type", 
                       tabName = "mv_type_tab"),
           menuSubItem("Choose Imputation Method", 
                       tabName = "imputation_tab")
  ),
  menuItem("Normalization", tabName = "nomalization",
           menuSubItem("Visualize Before Normalization",
                       tabName = "before_norm_tab"),
           menuSubItem("Choose Normalization Method ",
                        tabName = "normalization_method_tab"),
           menuSubItem("Visualize After Normalization",
                       tabName = "after_norm_tab")
  ),
  
        
  menuItem("Variation Inspection & Correction", tabName = "tmt_analysis",
           menuSubItem("Upload Metadata",
                       tabName = "metadata_upload_tab"),
           menuSubItem("Identify Sources of Variations",
                       tabName = "variation_source_tab"),
           menuSubItem("Batch Effect Inspection ",
                       tabName = "batch_inspection_tab"),
           menuSubItem("Batch Correction",
                       tabName = "batch_correction_tab"),
           menuSubItem("Visualize after batchcorrection",
                       tabName = "visualization_tab"),
           menuSubItem("check other variations",
                       tabName = "other_variations_tab")
  )
)

body <- dashboardBody(
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")),
  
  tabItems(
    tabItem(tabName = "tmt_data_upload_tab",
            fluidRow(
              box(width = 2,
                  title = "Upload TMT ProteinGroups",
                  solidHeader = T, status = "primary",
                  fileInput(inputId ="tmt_file", 
                            multiple = FALSE,
                            label = "",
                            accept = c("text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv"))
              ),
              
              box(width = 10, 
                  title = "Uploaded ProteinGroups File",
                  solidHeader = T, status = "primary",
                  dataTableOutput(outputId =  "tmt_raw_data")
              )
            )
    ),
    tabItem(tabName ="tmt_reportercorrected_tab",
            box(width = 2, 
                title = "Choose Working Variables",
                solidHeader = T, status = "primary",
                radioButtons(inputId = "reporter_intesi_corrected",
                             label = "",
                             choices = c("Reporter Intensity Corrected", "Reporter Intensity"))),
            box(width = 10, title = textOutput(outputId = "corrected_title"), 
                solidHeader = T, status = "primary",
                dataTableOutput(outputId =  "tmt_corrected_intensity"),
                downloadButton(outputId = "down_corrected_intensity_tmt")
            ),
           
            ),
    tabItem(tabName ="tmt_reporterintensity_tab",
            box(width = 9, title = "Filter Reporter Intensity Corrected", 
                solidHeader = T, status = "primary",
                dataTableOutput(outputId =  "tmt_reporter_intensity"),
                downloadButton(outputId = "down_reporter_intensity_tmt")
            )
            ),
    
    
    tabItem(tabName = "lfq_data_upload_tab",
            fluidRow(
              box(width = 2,
                  title = "Upload LFQ ProteinGroups",
                  solidHeader = T, status = "primary",
                  fileInput(inputId ="lfq_file", 
                            multiple = FALSE,
                            label = "",
                            accept = c("text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv"))
              ),
              
              box(width = 10, 
                  title = "Uploaded ProteinGroups File",
                  solidHeader = T, status = "primary",
                  dataTableOutput(outputId =  "lfq_raw_data")
              )
            )
    
  ),
  tabItem(tabName ="lfq_intensity_tab",
          box(width = 12, title = "LFQ Data Preparation", 
              solidHeader = T, status = "primary",
              dataTableOutput(outputId =  "lfq_intensity_data"),
              downloadButton(outputId = "down_lfq_intensity_lfq")
          )
  ),
  
 
  tabItem(tabName ="filter_missing_tab",
          
          box(width = 2, 
              title = "Adjust Threshold",
              solidHeader = T, status = "primary",
              sliderInput(inputId = "filter_threshold",
                           label = "",
                           min = 0,max = 100,value = 20)),
          box(width = 9, title = "Threshold Filtered Data", 
              solidHeader = T, status = "maroon",background = "aqua",
              dataTableOutput(outputId =  "threshold_data"),
              downloadButton(outputId = "down_threshold_data")
          )
  ),
  tabItem(tabName ="contaminants_filter_tab",
          box(width = 2, 
              title = "Choose contaminants",
              solidHeader = T, status = "primary",
              selectInput(inputId = "filter_contaminants",
                          label = " ",
                      choices = proteomics_contaminats,multiple = TRUE,
                      selected = "REV")),
          box(width = 10, title = "Filter Relevant Contaminants", 
              solidHeader = T, status = "maroon",background = "aqua",
              dataTableOutput(outputId =  "cont_free_data"),
              downloadButton(outputId = "down_cont_free_data")
          )
  ),
  tabItem(tabName ="mv_inspection_tab",
          box(width = 12, title = "Showing the Percentage of Missing Values in your Data", 
              solidHeader = T, status = "maroon",background = "aqua",
              dataTableOutput(outputId =  "missing_values_inspection",height = "500px"),
              downloadButton(outputId = "down_inspect_mv")
          )
  ),
  
  
  tabItem(tabName ="mv_type_tab",
          box(width = 12, title = "Type of Missing Values in  Your Data", 
              solidHeader = T, status = "maroon",
              htmlOutput(outputId =  "missing_values_type")
          )
  ),
  tabItem(tabName ="imputation_tab",
          box(width = 2, 
              title = "Choose Imputation Method",
              solidHeader = T,status = "maroon",
              radioButtons(inputId = "choose_imputation_method",
                           label = "",
                           choices = c("KNN", "missRanger","missForest","QRILC"))),
          box(width = 10, title = textOutput(outputId = "imputation_title"), 
              solidHeader = T, status = "maroon",background = "aqua",
              dataTableOutput(outputId =  "imputation_method"),
              downloadButton(outputId = "down_imputation_method")
          )
          
  ),
  tabItem(tabName ="before_norm_tab",
          box(width = 12, title = "Visualize the Data Distribution Before Normalization", 
              solidHeader = T, status = "maroon",background = "aqua",
              plotOutput(outputId =  "before_norm",height = "500px")
          )
  ),
  tabItem(tabName ="normalization_method_tab",
          box(width = 2, 
              title = "Choose Normalization Method",
              solidHeader = T, status = "maroon",
              radioButtons(inputId = "normalization_method",
                           label = "",
                           choices = c("VSN", "Log Normalization","Cyclicloess","Quantile"))),
        
          
          box(width = 10, title = textOutput(outputId ="normalization_title"),
              solidHeader = T, status = "maroon",background = "aqua",
              dataTableOutput(outputId =  "norm_method"),
              downloadButton(outputId = "down_norm_method")
          ),

          
  ),
  tabItem(tabName ="after_norm_tab",
          box(width = 12, title = "Visualize the Data Distribution After Normalization", 
              solidHeader = T, status = "maroon",background = "aqua",
            plotOutput(outputId =  "after_norm",height = "600px")
          )
  ),
  tabItem(tabName = "metadata_upload_tab",
          fluidRow(
            box(width = 2,
                title = "Upload Metadata file",
                solidHeader = T, status = "maroon",
                fileInput(inputId ="metadata_file", 
                          multiple = FALSE,
                          label = "",
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv"))
            ),
            
            box(width = 10,
                title = "Uploaded Metadata File",
                solidHeader = T, status = "maroon",
                dataTableOutput(outputId =  "metadata_data_df")
            )
          )
  ),
  tabItem(tabName ="variation_source_tab",
          box(width = 12, title = "Sources of Variations", 
              solidHeader = T, status = "maroon",
              plotOutput(outputId =  "variation_source_plot",height = "500px")
             
          )
  ),
  tabItem(tabName ="batch_inspection_tab",
          box(width = 12, title = "Confirm there is Batch Effect", 
              solidHeader = T, status = "maroon",
              plotOutput(outputId =  "inspect_batcheffect",height = "500px"),
              plotOutput(outputId = "down_inspect_batcheffect")
          )
  ),
  
  tabItem(tabName ="batch_correction_tab",
          box(width = 12, title = "Batch Correction using comBat", 
              solidHeader = T, status = "maroon",
              dataTableOutput(outputId =  "batch_correction"),
              downloadButton(outputId = "down_batch_correction")
          ),
          
  ),
  tabItem(tabName ="visualization_tab",
          box(width = 12, title = "Visualize After Batch Correction", 
              solidHeader = T, status = "maroon",
              plotOutput(outputId =  "visualize_after_correction",height = "500px")
          )
  ),
  tabItem(tabName ="other_variations_tab",
          box(width = 12, title = "Other Sources of Variations", 
              solidHeader = T, status = "maroon",
              plotOutput(outputId =  "other_variation_source",height = "500px"),
              plotOutput(outputId = "down_other_variation_source")
          )
  )

)
)            
            
            
# header$children[[2]]$children <-  tags$a(tags$img(src='ppw caption.PNG',height='100',width='600'))


ui <- dashboardPage(header,
                    sidebar,
                    body)

