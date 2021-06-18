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
require("dplyr") || utils::install.packages("dplyr")
require("stringr") || utils::install.packages("stringr")

#The following bioconductor packages are also needed for this application to work. If you haven't got them installed you can uncomment the following lines of code.

#BiocManager::install("Biostrings")
#BiocManager::install("GEOquery")
#BiocManager::install("SummarizedExperiment")

# Load R packages
library(tidyverse)
library(shiny)
library(shinythemes)
library(shinydashboard)
library(GEOquery)
library(DT)
library(pheatmap)
library(here)
library(dplyr)
library(SummarizedExperiment)
library(stringr)

#This file contains the server of the application.


shinyServer(function(input, output) {
    
    
    ## INSTRUCTIONS TAB server functionality now follows --------------------------------------------------------------------------------------------=--------------------------------------------------
    
    output$image1 <- renderImage({
        

        imagefile <- here::here("images", "Instructions_1.png")

        list(src = imagefile,
             width = 1000)
        
    }, deleteFile = FALSE)

    
    
    
    ## DATA DOWNLOADING TAB server functionality now follows    ----------------------------------------------------------------------------------------------------------------------------------------  
    gse_number <- eventReactive(input$action_gse_download, {
        input$txt_gse_download
    })
    
    #download the dataset
    gse_dataset <- reactive({
        gse_data_download_function(gse_number())
    }) 
    
    
    #output the contents of the gse dataset (shows different platforms)
    output$txtout_gse_download <- renderPrint(print(gse_dataset()))
    
    #phenodata of platform choice gets loaded into a variable
    phenodata <- eventReactive(input$action_select_platform, {
        gse_dataset()[[input$select_platform]] %>% 
            pData()
    })    
    
    
    output$tbl_platform_preview <- DT::renderDataTable(
        phenodata() %>%
            select(geo_accession, organism_ch1, characteristics_ch1, supplementary_file_1)
        ,
        style = "bootstrap",
        editable = TRUE,
        server = TRUE
    )
    
    
    ## SUPP DATA DOWNLOADING TAB server functionality now follows:    -----------------------------------------------------------------------------------------------------------------------------------
    

    supp_files <- eventReactive(input$action_download_supp_files, {
        gse_supp_download_function(gse_number())
        })
    
    output$list_of_downloaded_files <- renderTable(supp_files())
    
    
    observeEvent(input$action_untar_unzip_supp_files, {
        untar_unzip_function(gse_number())
    })
    
    output$list_of_untarred_files <- renderText(list.files(here::here("data", gse_number())))
    
    #render checkbox options sample choices
    output$sample_choices <- renderUI({
        checkboxGroupInput("selected_samples", label = NULL, choices = map(list.files(here::here("data", gse_number()), full.names = TRUE), basename))
    })
    
    ## GROUP SELECTION TAB server functionality now follows:  -------------------------------------------------------------------------------------------------------------------------------------------  
    
    #get a list of the column names present in the phenodata so the user can make a selection. 
    phenodata_columns <- reactive(colnames(phenodata()))
    
    #render checkbox options
    output$column_choices <- renderUI({
        checkboxGroupInput("selected_columns", label = NULL, choices = phenodata_columns())
    })
    
    #show columns of phenodata selected by userinput
    output$tbl_phenodata_preview <- DT::renderDataTable(
        phenodata() %>%
            dplyr::select(input$selected_columns),
        style = "bootstrap",
        server = TRUE
    )
    
    ## BUILD SUMMARIZED EXPERIMENT TAB server functionality now follows:   --------------------------------------------------------------------------------------------------------------------------------
    featuredata <- reactive({gse_dataset()[[input$select_platform]] %>% 
            featureData()
    })
    
    assaydata <- reactive({gse_dataset()[[input$select_platform]] %>% 
            exprs()
    })
    #Waarom staat deze code in het paars?    
    #gse <- reactive({
    #    SummarizedExperiment(assays = assaydata(),
    #                         rowData = featuredata()[[1]],
    #                         colData = phenodata())
    #})
    
    #output$tbl_assaydata_preview <- DT::renderDataTable(
    #    SummarizedExperiment::assay(se()),
    #    style = "bootstrap",
    #    server = TRUE
    #)
    
    se <- reactive({
        SummarizedExperiment(assays = assaydata(),
                             rowData = featuredata()[[1]],
                             colData = phenodata())
    })
    
    output$tbl_assaydata_preview <- DT::renderDataTable({
        input$build_sum

        (SummarizedExperiment::assay(se()))
         

    })
    
    ## testing tab server functionality now follows ------------------------------------------------------------------------------------------------------------------------------------------------------------
    ntext <- eventReactive(input$goButton, {
        input$n
    })
    
    output$nText <- renderText({
        ntext()
    })
    
    
    
})
