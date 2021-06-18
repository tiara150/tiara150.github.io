# Welcome to the differential expression analysis dashboard. 
# This application aims to make it easier for scientists to analyse RNA-seq data using the DESeq2 workflow developed by Michael I. Love, Simon Anders, and Wolfgang Huber.

require("tidyverse") || utils::install.packages("tidyverse")
require("shiny") || utils::install.packages("shiny")
require("shinythemes") || utils::install.packages("shinythemes")
require("GEOquery") || utils::install.packages("GEOquery")
require("DT") || utils::install.packages("DT")
require("pheatmap") || utils::install.packages("pheatmap")
require("here") || utils::install.packages("here")
require("shinydashboard") || utils::install.packages("shinydashboard")
require("stringr") || utils::install.packages("stringr")

#The following bioconductor packages are also needed for this application to work. If you haven't got them installed you can uncomment the following lines of code.

#BiocManager::install("Biostrings")
#BiocManager::install("GEOquery")
#BiocManager::install("SummarizedExperiment")

# Load R packages
library(tidyverse)
library(shiny)
library(shinydashboard)
library(shinythemes)
library(GEOquery)
library(DT)
library(pheatmap)
library(here)
library(SummarizedExperiment)
library(stringr)

#This file contains the userinterface of the application.

dashboardPage(      
  
  dashboardHeader(title = "DESeq2 Analysis"),
  
  dashboardSidebar(
    #make the menu -------------------------------------------------------------------------------------------------------------------------------------------------
    sidebarMenu(
      menuItem("Instructions", tabName = "instructions", icon = icon("info")),
      menuItem("Download Dataset", tabName = "download_dataset", icon = icon("download", lib = "font-awesome")),
      menuItem("Download Supp Files", tabName = "supp_files_assay_data", icon = icon("download", lib = "font-awesome")),
      menuItem("Select Experimental Groups", tabName = "select_groups", icon = icon("check", lib = "font-awesome")),
      menuItem("Build Summarized Experiment", tabName = "build_summarized_experiment", icon = icon("wrench", lib = "font-awesome")),
      menuItem("Testing Tab", tabName = "test_tab", icon = icon("vial", lib = "font-awesome"))
      )),
  
  
  dashboardBody(
    tabItems(
      #instructions tab --------------------------------------------------------------------------------------------------------------------------------------------
      tabItem(tabName = "instructions",
              
              titlePanel(tags$h1("Welcome to the DESeq2 Analysis Dashboard!")),
              
              fluidRow(
                box(
                  "instructions and tutorials go here."
                )),
              
              absolutePanel(
                imageOutput("image1")
              )
              
              
      ), 
      
      # data downloading tab ----------------------------------------------------------------------------------------------------------------------------------------
      tabItem(tabName = "download_dataset",
              fluidRow(
                #inputbox for the GSE number to download the dataset.
                box(
                  title = "GEO dataset Downloading",
                  width = "2",
                  textInput("txt_gse_download", "Enter GEO accession number. macaca dataset: GSE152439, datasets met featuredata: GSE153750, GSE142018, GSE140073"),
                  actionButton("action_gse_download", "Download this GSE dataset", icon("download", lib = "font-awesome"))
                ),
                
                #outputbox showing the platforms that are in the dataset
                box(
                  title = "These platforms are in the dataset",
                  width = "10",
                  verbatimTextOutput("txtout_gse_download")
                ),
                
                #inputbox to select the platform 
                box(
                  title = "Select the platform you want to use for this analysis:",
                  width = "2",
                  numericInput("select_platform", label = NULL, 1, 10, 1),
                  actionButton("action_select_platform","load this (sub)dataset", icon("refresh"))
                ),
                
                #outputbox showing just a bit of the phenodata of the platform so user can decide what platform. 
                box(
                  title = "Small preview of the data that is inside the platform",
                  width = "10",
                  dataTableOutput("tbl_platform_preview")
                )
                
              )),
      
      
      #  supplementary files download + make assaydata ------------------------------------------------------------------------------------------------------------------
      tabItem(tabName = "supp_files_assay_data",
              fluidRow(
                
                #select what columns to view
                box(
                  title = "download the supplementary files",
                  width = "10",
                  actionButton("action_download_supp_files", label = "download supp files", icon("download", lib = "font-awesome"))
              ),
              
              #show a list of supp files found in the folder
              box(
                title = "Show a list of supplementary files that have been downloaded.",
                width = "10",
                tableOutput("list_of_downloaded_files")
              ),
              
              #untar and unzip 
              box(
                title = "untar and unzip files",
                width = "10",
                actionButton("action_untar_unzip_supp_files", label = "untar and unzip supp files", icon("download", lib = "font-awesome"))
              ),
              
              #show files that have been untarred and unzipped
              box(
                title = "Show a list of supplementary files that have been unzipped.",
                width = "10",
                verbatimTextOutput("list_of_untarred_files")
              ),
              
              #sample selection box. uses names of files that are in: data --> "GSE###### folder" might need an actionbutton later depending on how the next steps are build. the boxes that are selected will be saved to input$sample_choices
              box(
                title = "Select samples",
                width = "3",
                uiOutput("sample_choices")
              ),
              
              fluidRow(
                box(
                  width = 10,
                  "Is the sample selection box empty? This is likely due to the folder being empty at the time the GSE number was entered. To fix this simply go to the 'download dataset' tab and click another time on the 'download this GSE dataset' button. 
                  The sample selection boxes should now appear. To fix this problem some sort of 'refresh' button needs to be coded but this is not very straight-forward when submitButtons can not be used."
                ))
              
        
             )),
      
    
      
      
      #  group selection tab ---------------------------------------------------------------------------------------------------------------------------------------
      # This tab used give users an option to select columns and show different parts of the coldata. Did not work anymore after rewriting the whole app to only use actionButtons. May or may not still be useful later.
      tabItem(tabName = "select_groups",
              fluidRow(
              
                #select what columns to view
                #box(
                #  title = "Select what columns of the phenodata you would like to view. This does not change the dataset but is just here for your convenience.",
                #  width = "2",
                #  uiOutput("column_choices")
                #),
                
                #show selected cols of phenodata
                #box(
                #  title = "Preview of the phenodata of this experiment",
                #  width = "10",
                #  dataTableOutput("tbl_phenodata_preview"),
                #  submitButton("Show selected columns", icon("refresh"))
                #)
                
              )),
      
      
      
      #  Build Summarized Experiment tab ---------------------------------------------------------------------------------------------------------------------------
      #this tab works when a dataset contains assaydata and makes a summarized experiment. Then shows a preview of the assaydata in the summarized experiment. 
      
    
      tabItem(tabName = "build_summarized_experiment",
              verticalLayout(
                
                #build experiment
                
                actionButton("build_sum", "build se", icon = icon("wrench", lib = "font-awesome")),
                
                
                #show selected cols of phenodata
                box(
                  title = "Preview of the assaydata of this summarized experiment",
                  width = "12",
                  dataTableOutput("tbl_assaydata_preview")
                )
                
              )),
      
      #  testing functions tab----------------------------------------------------------------------------------------------------------------------------------------
      # This is just a tab for the developer to test out functions and buttons etc before implementation in the app in the correct tab
      
      tabItem(tabName = "test_tab",
              fluidRow(
                pageWithSidebar(
                  headerPanel("actionButton test"),
                  sidebarPanel(
                    numericInput("n", "N:", min = 0, max = 100, value = 50),
                    br(),
                    actionButton("goButton", "Go!"),
                    p("Click the button to update the value displayed in the main panel.")
                  ),
                  mainPanel(
                    verbatimTextOutput("nText")
                  )
                )
           
             
              
              )
             )
      )
      
      
  )

)

