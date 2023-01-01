library(shiny)
library(tidyverse)
library(data.table)
library(ggpubr)
library(bslib)

light <- bs_theme(bootswatch = "minty")
dark <- bslib::bs_theme(
  bg = "#101010",
  fg = "#FDF7F7",
  primary = "#ED79F9",
  secondary = "#1e90ff"
)

ui <- fluidPage(
  theme = light,
  titlePanel(
    h1("Topological trajectory between clusters"),
    windowTitle = "Topological Trajectory - Shiny App"
  ),
  p(""),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      p("Explore the most likely topological reconfiguration trajectory between two age clusters"),
      checkboxInput("dark_mode", "Dark mode", value = FALSE),
      p("Visit", em(a("github.com/ClementGuichet/RestingState_LANG_Connectome", href = "https://github.com/ClementGuichet/RestingState_LANG_Connectome")), "for more information"),
      checkboxGroupInput("RSN_choice",
        h4("Pick resting-state networks"),
        choices = list(
          "Auditory" = "Auditory",
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
      checkboxInput("all", "Plot all regions"),
      selectInput("functional_role",
        h4("Explore modular or interareal functional role reconfiguration"),
        choices = list(
          "Modular" = "Modular",
          "Interareal" = "Interareal"
        )
      ),
      actionButton("button", "Plot trajectory", width = 460, style = "font-size: 20px"),
    ),
    mainPanel(
      h3(textOutput("select_var"), align = "center"),
      shinycssloaders::withSpinner(
        plotOutput("plot", height = 650),
        hide.ui = FALSE
      )
    )
  )
)


# Define server logic ----
server <- function(input, output, session) {
  observe(session$setCurrentTheme(
    if (isTRUE(input$dark_mode)) dark else light
  ))

  observe({
    x <- input$all

    if (!(is.null(x))) {
      y <- character(1)

      updateCheckboxGroupInput(session, "RSN_choice",
        choices = list(
          "Auditory" = "Auditory",
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
        selected = y
      )
    }
  })

  observe({
    x <- input$RSN_choice

    if (!(is.null(x))) {
      y <- character(1)

      updateCheckboxInput(session, "all", "Plot all regions", value = y)
    }
  })

  plot_button <- eventReactive(input$button, {
    input$RSN_choice
  })

  plot_button_bis <- eventReactive(input$button, {
    input$functional_role
  })

  output$select_var <- renderText({
    if (length(plot_button()) > 1) {
      list_input <- list()
      for (i in 1:length(plot_button())) {
        list_input[[i]] <- plot_button()[i]
        tmp_bis <- rbindlist(lapply(list_input, as.data.table))
        RSN_text <- capture.output(cat(tmp_bis %>% as.matrix(), sep = "|"))
      }
      paste("This is the most probable trajectory for", RSN_text, "regions combined.")
    } else if (length(plot_button()) == 1) {
      list_input <- list()
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
    modular <- read.csv("modular_trajectory.csv") %>%
      dplyr::select(-X) %>%
      plyr::rename(c("X1st_network" = "1st_network"))
    interareal <- read.csv("interareal_trajectory.csv") %>%
      dplyr::select(-X) %>%
      plyr::rename(c("X1st_network" = "1st_network"))

    trajectory <- function(type_func_df, cluster1, cluster2, list_RSN) {
      # Convert an adjacency dataframe to a 2-column dataframe
      adjacency_to_2col <- function(dataframe) {
        crossdata <- lapply(rownames(data), function(x) sapply(colnames(data), function(y) list(x, y, data[x, y])))
        crossdatatmp <- matrix(unlist(crossdata), nrow = 3)
        crossdatamat <- t(crossdatatmp)
        colnames(crossdatamat) <- c(cluster1, cluster2, "Value")
        crossdatadf <- as.data.frame(crossdatamat, stringsAsFactors = F)
        crossdatadf[, 3] <- as.numeric(crossdatadf[, 3])
        return(crossdatadf)
      }
      
      # For each region, find the most probable trajectory
      if (colnames(type_func_df[4]) == "Hub_consensus") {
        Cond_PMF <- type_func_df %>%
          spread(Age_group, freq) %>%
          arrange(Hub_consensus) %>%
          group_by(Region, .add = TRUE) %>%
          group_split()
      } else {
        Cond_PMF <- type_func_df %>%
          spread(Age_group, freq) %>%
          arrange(Bridgeness) %>%
          group_by(Region, .add = TRUE) %>%
          group_split()
      }
      
      
      Cond_PMF_list <- list()
      for (i in 1:length(Cond_PMF)) {
        tmp <- rbindlist(Cond_PMF[i])
        # Beware of the order here, first young then old
        data <- outer(tmp[, 5] %>% as.matrix(), tmp[, 4] %>% as.matrix()) %>% as.data.frame()
        if (colnames(type_func_df[4]) == "Hub_consensus") {
          colnames(data) <- c("Connector", "Peripheral", "Provincial", "Satellite")
          rownames(data) <- c("Connector", "Peripheral", "Provincial", "Satellite")
        } else {
          colnames(data) <- c("Global_Bridge", "Local_Bridge", "Not_a_Bridge", "Super_Bridge")
          rownames(data) <- c("Global_Bridge", "Local_Bridge", "Not_a_Bridge", "Super_Bridge")
        }
        crossdatadf <- adjacency_to_2col(data)
        data_2 <- crossdatadf %>% dplyr::slice_max(Value)
        Cond_PMF_list[[i]] <- cbind(tmp %>% slice(1:nrow(data_2)) %>% dplyr::select(`1st_network`, Region), data_2)
      }
      
      Cond_PMF_final <<- rbindlist(Cond_PMF_list) %>%
        # Identify each region with a unique label
        mutate(helper_vector = rep(seq(nrow(.)))) %>%
        unite(., Region, c("Region", "helper_vector"), remove = FALSE) %>%
        # Transform  back to an alluvial format
        pivot_longer(
          cols = c(cluster1, cluster2),
          names_to = "Age_group",
          values_to = colnames(type_func_df[4])
        )
      
      library(ggalluvial)
      if (list_RSN != "All") {
        if (colnames(type_func_df[4]) == "Hub_consensus") {
          display_cluster <- Cond_PMF_final %>%
            filter(grepl(list_RSN, `1st_network`)) %>%
            group_by(Age_group, Hub_consensus) %>%
            summarize(s = n()) %>%
            arrange(desc(Age_group), desc(Hub_consensus)) %>%
            .$Hub_consensus
          
          alluvial_cluster <- ggplot(
            Cond_PMF_final %>%
              filter(grepl(list_RSN, `1st_network`)),
            aes(x = forcats::fct_rev(Age_group), stratum = Hub_consensus, alluvium = Region, fill = Hub_consensus)
          ) +
            geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
            geom_stratum(alpha = .8) +
            scale_x_discrete(expand = c(.1, .1)) +
            # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
            geom_text(stat = "stratum", label = display_cluster) +
            scale_fill_brewer(palette = "Oranges", direction = 1) +
            labs(x = "Age group") +
            theme_pubclean()
          
          alluvial_cluster
        } else {
          display_cluster <- Cond_PMF_final %>%
            filter(grepl(list_RSN, `1st_network`)) %>%
            group_by(Age_group, Bridgeness) %>%
            summarize(s = n()) %>%
            arrange(desc(Age_group), desc(Bridgeness)) %>%
            .$Bridgeness
          
          alluvial_cluster <- ggplot(
            Cond_PMF_final %>%
              filter(grepl(list_RSN, `1st_network`)),
            aes(x = forcats::fct_rev(Age_group), stratum = Bridgeness, alluvium = Region, fill = Bridgeness)
          ) +
            geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
            geom_stratum(alpha = .8) +
            scale_x_discrete(expand = c(.1, .1)) +
            # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
            geom_text(stat = "stratum", label = display_cluster) +
            scale_fill_brewer(palette = "Oranges", direction = 1) +
            labs(x = "Age group") +
            theme_pubclean()
          
          alluvial_cluster
        }
      } else {
        if (colnames(type_func_df[4]) == "Hub_consensus") {
          display_cluster <- Cond_PMF_final %>%
            group_by(Age_group, Hub_consensus) %>%
            summarize(s = n()) %>%
            arrange(desc(Age_group), desc(Hub_consensus)) %>%
            .$Hub_consensus
          
          alluvial_cluster <- ggplot(
            Cond_PMF_final,
            aes(x = forcats::fct_rev(Age_group), stratum = Hub_consensus, alluvium = Region, fill = Hub_consensus)
          ) +
            geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
            geom_stratum(alpha = .8) +
            scale_x_discrete(expand = c(.1, .1)) +
            # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
            geom_text(stat = "stratum", label = display_cluster) +
            scale_fill_brewer(palette = "Oranges", direction = 1) +
            labs(x = "Age group") +
            theme_pubclean()
          
          alluvial_cluster
        } else {
          display_cluster <- Cond_PMF_final %>%
            group_by(Age_group, Bridgeness) %>%
            summarize(s = n()) %>%
            arrange(desc(Age_group), desc(Bridgeness)) %>%
            .$Bridgeness
          
          alluvial_cluster <- ggplot(
            Cond_PMF_final,
            aes(x = forcats::fct_rev(Age_group), stratum = Bridgeness, alluvium = Region, fill = Bridgeness)
          ) +
            geom_flow(alpha = .7, curve_type = "arctangent", width = .2, na.rm = TRUE) +
            geom_stratum(alpha = .8) +
            scale_x_discrete(expand = c(.1, .1)) +
            # geom_text(stat = "stratum", label = display_percentage, nudge_x = -0.05) +
            geom_text(stat = "stratum", label = display_cluster) +
            scale_fill_brewer(palette = "Oranges", direction = 1) +
            labs(x = "Age group") +
            theme_pubclean()
          
          alluvial_cluster
        }
      }
    }
    
    data <- switch(plot_button_bis(),
      "Modular" = modular,
      "Interareal" = interareal
    )

    if (length(plot_button()) == 1) {
      RSN <- input$RSN_choice
    } else if (length(plot_button()) > 1) {
      list_input <- list()
      for (i in 1:length(plot_button())) {
        list_input[[i]] <- plot_button()[i]
        tmp_bis <- rbindlist(lapply(list_input, as.data.table))
        RSN <- capture.output(cat(tmp_bis %>% as.matrix(), sep = "|"))
      }
    } else {
      RSN <- "All"
    }
    trajectory(data, "Young", "Old", RSN)
  })
}

shinyApp(ui = ui, server = server)
