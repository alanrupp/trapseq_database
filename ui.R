# Putting TRAP-seq data into an app that is organized and searchable
library(shiny)
library(dplyr)
library(DT)
library(shinythemes)


# - UI ------------------------------------------------------------------------
ui <- navbarPage(title = "TRAP-seq data", 
                 theme = shinytheme("flatly"),
                 
  # - Tab 1: Grab a dataset and view ------------------------------------------
  tabPanel("Explore data",
    fluidPage(
      
      # HTML tag info
      tags$head(
        tags$style(HTML("hr {border-top: 1px solid #000000;}"))
      ),
      
      # select dataset
      fluidRow(
        column(width = 3, 
          selectInput("pulldown", "Pulldown:", choices = pulldown_list)),
        column(width = 6,
               selectizeInput("biotypes", "Filter by biotype:",
                              choices = biotypes_list, multiple = TRUE,
                              selected = "protein_coding",
                              width = "100%")),
        column(width = 3, align = 'left', style = 'margin-top: 25px;',
               actionButton("select_data", "Select", width = "100%"))
        ),
      
      fluidRow(
        # plot all genes
        column(width = 6,
          plotOutput("plot1", width = "auto", height = "400px",
                     hover = hoverOpts("hover_data"))
        ),
        # output table of hovered genes
        column(width = 3,
               tableOutput("plot1_cursor"))
        ),
      
      # horizontal rule between plot and table for visual separation
      fluidRow(
        hr()
        ),
                                    
      # apply filters to the data that's displayed
      fluidRow(
        column(width = 3,
          checkboxInput("only_sig", "Significant (P < 0.05)", value = TRUE)),
        column(width = 3, 
          checkboxInput("only_enrich", "Enriched (>1.5 bead/sup)",
                        value = TRUE))
        ),
                                    
      # Display data table
      fluidRow(
        column(width = 9,
               DT::dataTableOutput("table")
        )
      ),
      
      # horizontal rule between table  and username for visual separation
      fluidRow(
        hr()
      ),
      tags$footer(
        textOutput("user")
      )
    )
    ),
  
  
  # - Tab 2: intersect datasets -----------------------------------------------
  tabPanel("Find unique genes",
    # Select input data
      fluidPage(
        sidebarLayout(
          sidebarPanel(
            selectizeInput("display_data", "Genes enriched in:", 
                           choices = pulldown_list, selected = pulldown_list[1]),
            selectizeInput("OR_data", "OR",
                           choices = pulldown_list,
                           multiple = TRUE),
            selectizeInput("AND_data", "AND",
                           choices = pulldown_list,
                           multiple = TRUE),
            selectizeInput("NOT_data", "NOT", 
                           choices = pulldown_list, 
                           multiple = TRUE),
            tags$hr(),
            numericInput("enrichment_threshold_b", "Enrichment threshold",
                         value = 1.5, min = 1),
            numericInput("bead_threshold_b", "Bead threshold (CPM)",
                         value = 1, min = 1),
            selectizeInput("biotypes_b", "Filter by biotype:",
                           choices = biotypes_list, multiple = TRUE),
            checkboxInput("dataset_search_sig", "Significant (P < 0.05)",
                          value = TRUE),
            actionButton("select_b", "Display", width = "100%")
          ),
          
          # display filtered genes
          mainPanel(DT::dataTableOutput("intersect_table"))
        )
      )
  ),
                 
  # - Tab 3: search by gene ---------------------------------------------------
  tabPanel("Find datasets by gene",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          selectizeInput("gene_search", "Gene", choices = gene_list, 
                         selected = "Lepr", multiple = FALSE),
          numericInput("gene_search_enrich", "Enrichment threshold", 
                       value = 1.5, min = 1),
          actionButton("select_gene", "Select"),
          checkboxInput("gene_search_sig", "Only Significant (P < 0.05)", 
                        value = TRUE)
        ),
       
      # List overlapping genes
        mainPanel(
          DT::dataTableOutput("datasets")
        )
      )
    )
  )
)