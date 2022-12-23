library(shiny)
library(tidyverse)
library(data.table)
library(ggpubr)

ui <- fluidPage(
  titlePanel("Topological trajectory between clusters", windowTitle = "Topological Trajectory - Shiny App"),
  
  sidebarLayout(
    sidebarPanel(
      p("The goal of this app is to display the most likely topological reconfiguration trajectory between two age clusters"),
      p("Visit", em(a("github.com/ClementGuichet/RestingState_LANG_Connectome", href="https://github.com/ClementGuichet/RestingState_LANG_Connectome")), "for more information"),
      checkboxGroupInput("RSN_choice", 
                         h4("Choose the resting-state networks to be displayed"),
        choices = list("Auditory" = "Auditory",
                       "CON" = "CON",
                       "DAN" = "DAN",
                       "DMN" = "DMN",
                       "FPN" = "FPN",
                       "Language" = "Language",
                       "SMN" = "SMN",
                       "PMM" = "PMM",
                       "VMM" = "VMM",
                       "Visual_1" = "Visual_1",
                       "Visual_2" = "Visual_2"
                       )
      ),
      p("You can freely choose to either visualize the", strong("modular or interareal"), "functional role reconfiguration"),
      selectInput("functional_role", 
                   h4("Choose the functional roles to be displayed"),
                   choices = list("Modular" = "Modular",
                                  "Interareal" = "Interareal"
                                  )),
      checkboxInput("all", "Plot all regions"),
      actionButton("button", "Plot trajectory", width = 150),
    ),
    mainPanel(
      h3(textOutput("select_var")),      
      shinycssloaders::withSpinner(
        plotOutput("plot", height = 650),
        hide.ui = FALSE
      )
    )
  )
)

# Define server logic ----
server <- function(input, output, session) {
  
  observe({
    x <- input$all
    
    if (is.null(x))
      x <- character(0)
    
    updateCheckboxGroupInput(session, "RSN_choice",
                             choices = list("Auditory" = "Auditory",
                                            "CON" = "CON",
                                            "DAN" = "DAN",
                                            "DMN" = "DMN",
                                            "FPN" = "FPN",
                                            "Language" = "Language",
                                            "SMN" = "SMN",
                                            "PMM" = "PMM",
                                            "VMM" = "VMM",
                                            "Visual_1" = "Visual_1",
                                            "Visual_2" = "Visual_2"
                             ),
                             selected = x
    )
  })
  
  plot_button <- eventReactive(input$button, {
    input$RSN_choice
  })
  
  plot_button_bis <- eventReactive(input$button, {
    input$functional_role
  })
  
  output$select_var <- renderText({
    if (length(plot_button()) > 1) {
      list_input <-  list()
      for (i in 1:length(plot_button())) {
        list_input[[i]] <- plot_button()[i]
        tmp_bis <- rbindlist(lapply(list_input, as.data.table))
        RSN_text <- capture.output(cat(tmp_bis %>% as.matrix(), sep = "|"))
      }
      paste("This is the most probable trajectory for", RSN_text, "regions combined.")
    } else if (length(plot_button()) == 1){
      list_input <-  list()
      for (i in 1:length(plot_button())) {
        list_input[[i]] <- plot_button()[i]
        tmp_bis <- rbindlist(lapply(list_input, as.data.table))
        RSN_text <- capture.output(cat(tmp_bis %>% as.matrix(), sep = "|"))
      }
      paste("This is the most probable trajectory for", RSN_text, "regions.")
    } else {
      paste("This is the most probable trajectory for all regions combined.")
    }
  })
  
  output$plot <- renderPlot({ 
    modular <- read.csv("modular_trajector.csv") %>% dplyr::select(-X) %>% plyr::rename(c("X1st_network" = "1st_network"))
    interareal <- read.csv("interareal_trajector.csv") %>% dplyr::select(-X) %>% plyr::rename(c("X1st_network" = "1st_network"))
    
    trajectory <- function(type_func_df, list_RSN) {
      # Convert an adjacency dataframe to a 2-column dataframe
      adjacency_to_2col <- function(dataframe) {
        crossdata <- lapply(rownames(data), function(x) sapply(colnames(data), function(y) list(x, y, data[x, y])))
        crossdatatmp <- matrix(unlist(crossdata), nrow = 3)
        crossdatamat <- t(crossdatatmp)
        colnames(crossdatamat) <- c("23", "56", "Value")
        crossdatadf <- as.data.frame(crossdatamat, stringsAsFactors = F)
        crossdatadf[, 3] <- as.numeric(crossdatadf[, 3])
        return(crossdatadf)
      }
      
      # For each region, find the most probable trajectory
      if (colnames(type_func_df[4]) == "Hub_consensus") {
        Cond_PMF <- type_func_df %>%
          spread(cluster, freq) %>%
          arrange(Hub_consensus) %>% 
          group_by(Region, .add = TRUE) %>%
          group_split()
      } else {
        Cond_PMF <- type_func_df %>%
          spread(cluster, freq) %>%
          arrange(Bridgeness) %>% 
          group_by(Region, .add = TRUE) %>%
          group_split()
      }
      
      
      Cond_PMF_list <- list()
      for (i in 1:length(Cond_PMF)) {
        tmp <- rbindlist(Cond_PMF[i])
        data <- outer(tmp$`23`, tmp$`56`) %>% as.data.frame()
        if (colnames(type_func_df[4]) == "Hub_consensus") {
          colnames(data) <- c("Connector", "Isolate", "Peripheral", "Provincial", "Satellite")
          rownames(data) <- c("Connector", "Isolate", "Peripheral", "Provincial", "Satellite")
        } else {
          colnames(data) <- c("Global_Bridge", "Local_Bridge", "Not_a_Bridge", "Super_Bridge")
          rownames(data) <- c("Global_Bridge", "Local_Bridge", "Not_a_Bridge", "Super_Bridge")
        }
        crossdatadf <- adjacency_to_2col(data)
        data_2 <- crossdatadf %>% dplyr::slice_max(Value)
        data_3 <- cbind(tmp %>% slice(1:nrow(data_2)) %>% dplyr::select(`1st_network`, Region), data_2)
        Cond_PMF_list[[i]] <- data_3
      }
      
      Cond_PMF_final <- rbindlist(Cond_PMF_list) %>%
        # Identify each region with a unique label
        mutate(helper_vector = rep(seq(nrow(.)))) %>%
        unite(., Region, c("Region", "helper_vector"), remove = FALSE) %>%
        # Transform  back to an alluvial format
        pivot_longer(
          cols = c("23", "56"),
          names_to = "cluster",
          values_to = colnames(type_func_df[4])
        )
      
      library(ggalluvial)
      
      if (colnames(type_func_df[4]) == "Hub_consensus") {
        display_cluster <- Cond_PMF_final %>%
          filter(grepl(list_RSN, `1st_network`)) %>% 
          group_by(cluster, Hub_consensus) %>%
          summarize(s = n()) %>%
          arrange(cluster, desc(Hub_consensus)) %>%
          .$Hub_consensus
        
        alluvial_cluster <- ggplot(
          Cond_PMF_final %>% 
            filter(grepl(list_RSN, `1st_network`)),
          aes(x = cluster, stratum = Hub_consensus, alluvium = Region, fill = Hub_consensus)
        ) +
          geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
          geom_stratum(alpha = .8) +
          scale_x_discrete(expand = c(.1, .1)) +
          # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
          geom_text(stat = "stratum", label = display_cluster) +
          scale_fill_brewer(palette = "Oranges", direction = 1) +
          theme_pubclean()
        
        alluvial_cluster
        
      } else {
        display_cluster <- Cond_PMF_final %>%
          filter(grepl(list_RSN, `1st_network`)) %>% 
          group_by(cluster, Bridgeness) %>%
          summarize(s = n()) %>%
          arrange(cluster, desc(Bridgeness)) %>%
          .$Bridgeness
        
        alluvial_cluster <- ggplot(
          Cond_PMF_final %>% 
            filter(grepl(list_RSN, `1st_network`)),
          aes(x = cluster, stratum = Bridgeness, alluvium = Region, fill = Bridgeness)
        ) +
          geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
          geom_stratum(alpha = .8) +
          scale_x_discrete(expand = c(.1, .1)) +
          # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
          geom_text(stat = "stratum", label = display_cluster) +
          scale_fill_brewer(palette = "Oranges", direction = 1) +
          theme_pubclean()
        
        alluvial_cluster
      }
    }
    
    data <- switch(plot_button_bis(),
                   "Modular" = modular,
                   "Interareal" = interareal)
    
    if (length(plot_button()) == 1) {
      RSN <- input$RSN_choice
    } else if (length(plot_button()) > 1) {
      list_input <-  list()
      for (i in 1:length(plot_button())) {
        list_input[[i]] <- plot_button()[i]
        tmp_bis <- rbindlist(lapply(list_input, as.data.table))
        RSN <- capture.output(cat(tmp_bis %>% as.matrix(), sep = "|"))
      }
    } else {
      RSN <- "Auditory|CON|DAN|DMN|FPN|Language|SMN|PMM|VMM|Visual_1|Visual_2"
    }
    trajectory(data, RSN)
  })
}

shinyApp(ui = ui, server = server)





