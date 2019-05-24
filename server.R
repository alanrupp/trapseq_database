library(dplyr)
library(ggplot2)
library(shiny)
library(DT)
library(readr)

# - Server function ----------------------------------------------------------
server <- function(input, output, session) {
  session$onSessionEnded(stopApp)
  
  output$user <- renderText(
    paste0("signed in as: ", 
           ifelse(is.null(session$user), "anonymous", session$user))
  )
  
  # read in data
  data_a <- eventReactive(input$select_data, {
    
    # grab data
    data_pop = input$pulldown
    data_area = input$region
    dataset <- read_csv(paste0("data/", input$region, "_", input$pulldown, ".csv"))
    
    # make sure requested data exists
    validate(
      need(nrow(dataset) > 0,
           message = c("Error -- Dataset does not exist. Try:\n",
                       c(paste(input$pulldown, areas[input$pulldown == pops]),
                         "or",
                         paste(pops[input$region == areas], input$region)
                       )
           )
      )
    )
    
    # arrange by Enrichment value
    dataset <- dataset %>%
      arrange(desc(P < 0.05), desc(Enrichment))
    
    
    # filter data by biotypes
    if (length(input$biotypes) != 0) {
      biotypes_selected <- biotypes_list[biotypes_list %in% input$biotypes]
      
      dataset <- dataset %>%
        inner_join(., genes, by = "Gene") %>%
        filter(Biotype %in% biotypes_selected) %>%
        select(-Biotype)
    }
    
    # export data, along with name of request, to a list
    list_a <- list(dataset, data_pop, data_area)
    names(list_a) <- c("dataset", "pop", "area")
    
    return(list_a)
    
  })
  
  # - Tab 1: Display table of data from a given experiment --------------------
  output$plot1 <- renderPlot({
    
    ggplot(data_a()[["dataset"]], 
           aes(x = Bead, y = Enrichment, color = P < 0.05)) +
      geom_hline(aes(yintercept = 0)) +
      geom_point(alpha = 0.4) +
      scale_color_manual(values = c("black", "red")) +
      scale_x_continuous(trans = "log2") +
      labs(title = paste(data_a()[["pop"]], data_a()[["area"]]),
           subtitle = "Hover over a point for info",
           x = "Bead expression (FPM)",
           y = "Enrichment (log2)") +
      theme_classic() +
      theme(legend.position = c(0.85, 0.15))
    
  })
  
  # table of points hovered over
  output$plot1_cursor <- renderTable({

    nearPoints(data_a()[["dataset"]], input$hover_data)
    
  })
  
  # output data table of raw data
  output$table <- renderDataTable({
    
   df <- data_a()[["dataset"]]
    
    if (input$only_sig == TRUE) {
      df <- filter(df, P < 0.05)
    }
    if (input$only_enrich == TRUE) {
      df <- filter(df, Enrichment >= 1.5)
    }
    
    return(df)
    
  })
  
  # -- Tab 2: List of overlapping genes in 2 (or more?) datasets --------------
  dsets <- eventReactive(input$select_b, {
    
    grab <- union(input$display_data,
                  union(input$OR_data,
                        union(input$AND_data, input$NOT_data)
                  )
    )
    
    # grab datasets and filter by enrichment and bead threshold
    df <- all_data[grab] %>%
      purrr::map(~dplyr::filter(.x,
        Enrichment >= input$enrichment_threshold_b &
        Bead >= input$bead_threshold_b)
        ) %>%
      purrr::set_names(grab)
    
    # filter data by biotypes
    if (length(input$biotypes_b) != 0) {
      biotypes_selected <- biotypes_list[biotypes_list %in% input$biotypes]
      
      df <- df %>%
        inner_join(., genes, by = "Gene") %>%
        filter(Biotype %in% biotypes_selected) %>%
        select(-Biotype)
    }
    
    # if only significant datasets:
    if (input$dataset_search_sig == TRUE) {
      df <- df %>%
        purrr::map(~dplyr::filter(.x, P < 0.05))
    }
    
    display_df <- df[[input$display_data]]
    or_df <- df[[input$OR_data]]
    and_df <- df[[input$AND_data]]
    not_df <- df[[input$NOT_data]]
    
    return(display_df)
    
  })

  
  # Display table of intersecting genes
  output$intersect_table <- renderDataTable(
    DT::datatable({
      
      dsets()
      
      })
    )
  
  
  # -- Tab 3: Dataset overlap by Gene -----------------------------------------
  gene <- eventReactive(input$select_gene, {
    
    # ensure gene is in gene list
    validate(
      need(input$gene_search %in% gene_list, 
           message = paste("Error -- invalid gene name")
      )
    )
    
    # grab the rows from all data specified by gene input
    df <- 
      dplyr::bind_rows(
      purrr::map(all_data, 
                 ~dplyr::filter(.x, 
                                Gene == as.character(input$gene_search))
                 )
      )
    
    # add column for dataset of origin and filter by specified enrich threshold
    df <- df %>%
      dplyr::mutate(Dataset = datasets) %>%
      dplyr::filter(Enrichment >= input$gene_search_enrich) %>%
      select("Dataset", "Gene", "Bead", "Enrichment", "P")
    
    if (input$gene_search_sig == TRUE) {
      df <- df %>%
        dplyr::filter(P < 0.05)
    }
    
    # ensure at least 1 dataset has gene enriched
    validate(
      need(nrow(df) > 0, 
           if (input$gene_search_sig == TRUE) {
             message = paste("No datasets have", input$gene_search, 
                             "significantly enriched at", 
                             input$gene_search_enrich, "fold-change")
             
           } else {
             message = paste("No datasets have", input$gene_search, 
                             "enriched at", 
                             input$gene_search_enrich, "fold-change")
           }
      )
    )
    
    # return data table
    return(df)
    
  })
  
  # display datasets that feature that gene
  output$datasets <- renderDataTable(
    DT::datatable({
      
      gene()
      
      }, 
      options = list(searching = FALSE)))
  
}
