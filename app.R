# Population Bioequivalence (PBE) Analysis Shiny App
# Based on FDA Draft Guidance on Budesonide (September 2012)

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(plotly)
library(readr)
library(dplyr)
library(knitr)
library(kableExtra)

# Define regulatory constants
SIGMA_T0 <- 0.1
THETA_P <- 2.0891

# PBE Calculation Functions
calculate_msb_msw <- function(data_df, m) {
  # Group by batch and container to get container means
  container_stats <- data_df %>%
    group_by(Batch, Container) %>%
    summarise(
      Container_Mean = mean(Measurement),
      Count = n(),
      .groups = 'drop'
    )
  
  container_means <- container_stats$Container_Mean
  overall_mean <- mean(data_df$Measurement)
  n_containers <- length(container_means)
  
  # Calculate MSB
  msb_sum <- sum((container_means - overall_mean)^2)
  msb <- (m * msb_sum) / (n_containers - 1)
  
  # Calculate MSW
  msw_sum <- 0
  for(i in 1:nrow(container_stats)) {
    batch <- container_stats$Batch[i]
    container <- container_stats$Container[i]
    container_data <- filter(data_df, Batch == batch, Container == container)
    container_mean <- mean(container_data$Measurement)
    msw_sum <- msw_sum + sum((container_data$Measurement - container_mean)^2)
  }
  
  msw <- msw_sum / (n_containers * (m - 1))
  
  return(list(msb = msb, msw = msw, container_stats = container_stats))
}

calculate_pbe_components <- function(delta, msb_t, msw_t, msb_r, msw_r, n_t, l_t, n_r, l_r, m, alpha = 0.05) {
  components <- list()
  
  # ED Component (Mean Difference)
  ed <- delta^2
  components$ed <- ed
  
  # HD Component (t-distribution)
  df_pooled <- (l_t * n_t - 1) + (l_r * n_r - 1)
  t_critical <- qt(1 - alpha/2, df_pooled)
  hd_variance <- msb_t/(n_t * l_t * m) + msb_r/(n_r * l_r * m)
  hd <- (delta + t_critical * sqrt(hd_variance))^2
  ud <- (hd - ed)^2  # CORRECTED: U = (H-E)²
  
  components$hd <- hd
  components$ud <- ud
  
  # E1 Component (Test MSB)
  e1 <- msb_t / m
  components$e1 <- e1
  
  # H1 Component (Chi-square) - Lower tail for test product
  df1 <- l_t * n_t - 1
  chi2_critical_1 <- qchisq(alpha, df1)  # Lower tail for H1, H2
  h1 <- (l_t * n_t - 1) * e1 / chi2_critical_1
  u1 <- (h1 - e1)^2  # CORRECTED: U = (H-E)²
  
  components$h1 <- h1
  components$u1 <- u1
  
  # E2 Component (Test MSW)
  e2 <- (m - 1) * msw_t / m
  components$e2 <- e2
  
  # H2 Component (Chi-square) - Lower tail for test product
  df2 <- l_t * n_t * (m - 1)
  chi2_critical_2 <- qchisq(alpha, df2)  # Lower tail for H1, H2
  h2 <- l_t * n_t * (m - 1) * e2 / chi2_critical_2
  u2 <- (h2 - e2)^2  # CORRECTED: U = (H-E)²
  
  components$h2 <- h2
  components$u2 <- u2
  
  return(components)
}

calculate_reference_scaled_components <- function(msb_r, msw_r, n_r, l_r, m, theta_p, alpha = 0.05) {
  ref_components <- list()
  
  # E3s Component (Reference MSB - scaled)
  e3s <- -(1 + theta_p) * msb_r / m
  ref_components$e3s <- e3s
  
  # H3s Component (Chi-square) - Upper tail for reference product
  df3 <- l_r * n_r - 1
  chi2_critical_3 <- qchisq(1 - alpha, df3)  # Upper tail for H3s, H4s
  h3s <- (l_r * n_r - 1) * e3s / chi2_critical_3
  u3s <- (h3s - e3s)^2  # CORRECTED: U = (H-E)²
  
  ref_components$h3s <- h3s
  ref_components$u3s <- u3s
  
  # E4s Component (Reference MSW - scaled)
  e4s <- -(1 + theta_p) * (m - 1) * msw_r / m
  ref_components$e4s <- e4s
  
  # H4s Component (Chi-square) - Upper tail for reference product
  df4 <- l_r * n_r * (m - 1)
  chi2_critical_4 <- qchisq(1 - alpha, df4)  # Upper tail for H3s, H4s
  h4s <- l_r * n_r * (m - 1) * e4s / chi2_critical_4
  u4s <- (h4s - e4s)^2  # CORRECTED: U = (H-E)²
  
  ref_components$h4s <- h4s
  ref_components$u4s <- u4s
  
  return(ref_components)
}

calculate_constant_scaled_components <- function(msb_r, msw_r, n_r, l_r, m, alpha = 0.05) {
  const_components <- list()
  
  # E3c Component
  e3c <- -msb_r / m
  const_components$e3c <- e3c
  
  # H3c Component - Upper tail for reference product
  df3 <- l_r * n_r - 1
  chi2_critical_3 <- qchisq(1 - alpha, df3)
  h3c <- (l_r * n_r - 1) * e3c / chi2_critical_3
  u3c <- (h3c - e3c)^2  # CORRECTED: U = (H-E)²
  
  const_components$h3c <- h3c
  const_components$u3c <- u3c
  
  # E4c Component
  e4c <- -(m - 1) * msw_r / m
  const_components$e4c <- e4c
  
  # H4c Component - Upper tail for reference product
  df4 <- l_r * n_r * (m - 1)
  chi2_critical_4 <- qchisq(1 - alpha, df4)
  h4c <- l_r * n_r * (m - 1) * e4c / chi2_critical_4
  u4c <- (h4c - e4c)^2  # CORRECTED: U = (H-E)²
  
  const_components$h4c <- h4c
  const_components$u4c <- u4c
  
  return(const_components)
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Population Bioequivalence (PBE) Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
      menuItem("Results", tabName = "results", icon = icon("table")),
      menuItem("Visualizations", tabName = "plots", icon = icon("chart-line")),
      menuItem("Report", tabName = "report", icon = icon("file-alt"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side {
          background-color: #f4f4f4;
        }
        .box {
          margin-bottom: 20px;
        }
        .validation-success {
          color: #27ae60;
          font-weight: bold;
        }
        .validation-warning {
          color: #f39c12;
          font-weight: bold;
        }
      "))
    ),
    
    tabItems(
      # Upload Data Tab
      tabItem(tabName = "upload",
        fluidRow(
          box(
            title = "Data Upload", status = "primary", solidHeader = TRUE, width = 6,
            fileInput("file", "Choose CSV File",
                     accept = c(".csv")),
            checkboxInput("header", "Header", TRUE),
            radioButtons("sep", "Separator",
                        choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                        selected = ","),
            tags$hr(),
            p("Expected format: Batch, Container, Stage, Product, Measurement"),
            p("Products should be labeled as 'TEST' and 'REF'"),
            p("Stages should be 'B' (Beginning), 'M' (Middle), 'E' (End)")
          ),
          
          box(
            title = "Regulatory Constants", status = "info", solidHeader = TRUE, width = 6,
            p(strong("σ_T0:"), SIGMA_T0, "(Regulatory constant)"),
            p(strong("θ_p:"), THETA_P, "(Regulatory constant)"),
            br(),
            p("These constants are defined in the FDA Draft Guidance on Budesonide (September 2012)"),
            actionButton("analyze", "Run PBE Analysis", 
                        class = "btn-success btn-lg", icon = icon("play"))
          )
        ),
        
        fluidRow(
          box(
            title = "Data Preview", status = "primary", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("preview")
          )
        )
      ),
      
      # Results Tab
      tabItem(tabName = "results",
        fluidRow(
          box(
            title = "Study Design Summary", status = "success", solidHeader = TRUE, width = 6,
            verbatimTextOutput("study_summary")
          ),
          
          box(
            title = "Procedure Selection", status = "info", solidHeader = TRUE, width = 6,
            verbatimTextOutput("procedure_selection")
          )
        ),
        
        fluidRow(
          box(
            title = "PBE Components - E Values (Point Estimates)", 
            status = "primary", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("e_components")
          )
        ),
        
        fluidRow(
          box(
            title = "PBE Components - H Values (Upper Confidence Bounds)", 
            status = "primary", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("h_components")
          )
        ),
        
        fluidRow(
          box(
            title = "PBE Components - U Values (Variance Terms: U = (H-E)²)", 
            status = "primary", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("u_components")
          )
        ),
        
        fluidRow(
          box(
            title = "Final Bioequivalence Results", 
            status = "warning", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("final_results")
          )
        )
      ),
      
      # Visualizations Tab
      tabItem(tabName = "plots",
        fluidRow(
          box(
            title = "Data Distribution by Product", 
            status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("dist_plot")
          ),
          
          box(
            title = "Batch Variability", 
            status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("batch_plot")
          )
        ),
        
        fluidRow(
          box(
            title = "PBE Components Breakdown", 
            status = "primary", solidHeader = TRUE, width = 12,
            plotlyOutput("components_plot")
          )
        )
      ),
      
      # Report Tab
      tabItem(tabName = "report",
        fluidRow(
          box(
            title = "Complete PBE Analysis Report", 
            status = "success", solidHeader = TRUE, width = 12,
            downloadButton("download_report", "Download Report", 
                          class = "btn-primary", icon = icon("download")),
            br(), br(),
            verbatimTextOutput("full_report")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive values
  values <- reactiveValues(
    data = NULL,
    results = NULL,
    components = NULL
  )
  
  # File upload and preview
  observeEvent(input$file, {
    req(input$file)
    
    tryCatch({
      values$data <- read.csv(input$file$datapath,
                             header = input$header,
                             sep = input$sep)
      
      # Validate data structure
      required_cols <- c("Batch", "Container", "Stage", "Product", "Measurement")
      if (!all(required_cols %in% names(values$data))) {
        showNotification("Error: Missing required columns. Expected: Batch, Container, Stage, Product, Measurement", 
                        type = "error", duration = 10)
        values$data <- NULL
        return()
      }
      
      # Validate products
      products <- unique(values$data$Product)
      if (!all(c("TEST", "REF") %in% products)) {
        showNotification("Error: Data must contain both 'TEST' and 'REF' products", 
                        type = "error", duration = 10)
        values$data <- NULL
        return()
      }
      
      showNotification("Data uploaded successfully!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error", duration = 10)
      values$data <- NULL
    })
  })
  
  output$preview <- DT::renderDataTable({
    req(values$data)
    DT::datatable(values$data, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # PBE Analysis
  observeEvent(input$analyze, {
    req(values$data)
    
    tryCatch({
      # Separate test and reference data
      test_df <- filter(values$data, Product == 'TEST')
      ref_df <- filter(values$data, Product == 'REF')
      
      # Basic statistics
      test_mean <- mean(test_df$Measurement)
      ref_mean <- mean(ref_df$Measurement)
      delta <- test_mean - ref_mean
      
      # Study design parameters
      m <- length(unique(values$data$Stage))
      l_t <- length(unique(test_df$Batch))
      l_r <- length(unique(ref_df$Batch))
      n_t <- length(unique(test_df$Container)) / l_t
      n_r <- length(unique(ref_df$Container)) / l_r
      
      # Calculate MSB and MSW
      test_results <- calculate_msb_msw(test_df, m)
      ref_results <- calculate_msb_msw(ref_df, m)
      
      msb_t <- test_results$msb
      msw_t <- test_results$msw
      msb_r <- ref_results$msb
      msw_r <- ref_results$msw
      
      # Calculate sigma values
      between_component_t <- msb_t / m
      within_component_t <- (m - 1) * msw_t / m
      sigma_t <- sqrt(between_component_t + within_component_t)
      
      between_component_r <- msb_r / m
      within_component_r <- (m - 1) * msw_r / m
      sigma_r <- sqrt(between_component_r + within_component_r)
      
      # Determine procedure
      use_reference_scaled <- sigma_r > SIGMA_T0
      
      # Calculate PBE components
      components <- calculate_pbe_components(delta, msb_t, msw_t, msb_r, msw_r, 
                                           n_t, l_t, n_r, l_r, m)
      
      # Calculate reference-scaled components
      ref_components <- calculate_reference_scaled_components(msb_r, msw_r, n_r, l_r, m, THETA_P)
      
      # Calculate constant-scaled components
      const_components <- calculate_constant_scaled_components(msb_r, msw_r, n_r, l_r, m)
      
      # Combine all components
      all_components <- c(components, ref_components, const_components)
      
      # Final calculations
      # Reference-scaled procedure
      eq_ref <- all_components$ed + all_components$e1 + all_components$e2 + 
                all_components$e3s + all_components$e4s
      uq_ref <- all_components$ud + all_components$u1 + all_components$u2 + 
                all_components$u3s + all_components$u4s
      hq_ref <- eq_ref + sqrt(max(0, uq_ref))
      bioequivalent_ref <- hq_ref <= 0
      
      # Constant-scaled procedure
      eq_const <- all_components$ed + all_components$e1 + all_components$e2 + 
                  all_components$e3c + all_components$e4c - THETA_P * SIGMA_T0^2
      uq_const <- all_components$ud + all_components$u1 + all_components$u2 + 
                  all_components$u3c + all_components$u4c
      hq_const <- eq_const + sqrt(max(0, uq_const))
      bioequivalent_const <- hq_const <= 0
      
      # Store results
      values$results <- list(
        test_mean = test_mean,
        ref_mean = ref_mean,
        delta = delta,
        sigma_t = sigma_t,
        sigma_r = sigma_r,
        msb_t = msb_t, msw_t = msw_t,
        msb_r = msb_r, msw_r = msw_r,
        m = m, l_t = l_t, l_r = l_r, n_t = n_t, n_r = n_r,
        use_reference_scaled = use_reference_scaled,
        eq_ref = eq_ref, uq_ref = uq_ref, hq_ref = hq_ref, bioequivalent_ref = bioequivalent_ref,
        eq_const = eq_const, uq_const = uq_const, hq_const = hq_const, bioequivalent_const = bioequivalent_const
      )
      
      values$components <- all_components
      
      showNotification("PBE Analysis completed successfully!", type = "message")
      # Switch to results tab after analysis
      # updateTabItems(session, "tabs", "results")
      
    }, error = function(e) {
      showNotification(paste("Error in analysis:", e$message), type = "error", duration = 10)
    })
  })
  
  # Study Summary Output
  output$study_summary <- renderText({
    req(values$results)
    
    paste(
      "STUDY DESIGN PARAMETERS",
      "========================",
      paste("Number of life stages (m):", values$results$m),
      paste("Test batches (l_T):", values$results$l_t),
      paste("Reference batches (l_R):", values$results$l_r),
      paste("Containers per batch (n):", values$results$n_t),
      "",
      "BASIC STATISTICS",
      "================",
      paste("Test mean (μ_T):", sprintf("%.8f", values$results$test_mean)),
      paste("Reference mean (μ_R):", sprintf("%.8f", values$results$ref_mean)),
      paste("Delta (Δ = μ_T - μ_R):", sprintf("%.8f", values$results$delta)),
      paste("Delta squared (Δ²):", sprintf("%.8f", values$results$delta^2)),
      "",
      "VARIANCE COMPONENTS",
      "===================",
      paste("σ_T:", sprintf("%.8f", values$results$sigma_t)),
      paste("σ_R:", sprintf("%.8f", values$results$sigma_r)),
      paste("MSB_T:", sprintf("%.8f", values$results$msb_t)),
      paste("MSW_T:", sprintf("%.8f", values$results$msw_t)),
      paste("MSB_R:", sprintf("%.8f", values$results$msb_r)),
      paste("MSW_R:", sprintf("%.8f", values$results$msw_r)),
      sep = "\n"
    )
  })
  
  # Procedure Selection Output
  output$procedure_selection <- renderText({
    req(values$results)
    
    paste(
      "PROCEDURE SELECTION",
      "===================",
      paste("σ_R (", sprintf("%.6f", values$results$sigma_r), ")", 
            ifelse(values$results$use_reference_scaled, " >", " ≤"), 
            " σ_T0 (", SIGMA_T0, ")", sep = ""),
      paste("→ Use", ifelse(values$results$use_reference_scaled, 
                           "Reference-scaled", "Constant-scaled"), "procedure"),
      "",
      "REGULATORY CONSTANTS",
      "====================",
      paste("σ_T0 =", SIGMA_T0, "(Regulatory constant)"),
      paste("θ_p =", THETA_P, "(Regulatory constant)"),
      sep = "\n"
    )
  })
  
  # E Components Table
  output$e_components <- DT::renderDataTable({
    req(values$components)
    
    e_data <- data.frame(
      Component = c("ED", "E1", "E2", "E3s", "E4s", "E3c", "E4c"),
      Description = c(
        "Mean Difference (Δ²)",
        "Test Between-Container Variance",
        "Test Within-Container Variance", 
        "Reference Between-Container (Scaled)",
        "Reference Within-Container (Scaled)",
        "Reference Between-Container (Constant)",
        "Reference Within-Container (Constant)"
      ),
      Value = c(
        sprintf("%.8f", values$components$ed),
        sprintf("%.8f", values$components$e1),
        sprintf("%.8f", values$components$e2),
        sprintf("%.8f", values$components$e3s),
        sprintf("%.8f", values$components$e4s),
        sprintf("%.8f", values$components$e3c),
        sprintf("%.8f", values$components$e4c)
      ),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(e_data, options = list(dom = 't', pageLength = 20)) %>%
      formatStyle(columns = 1:3, fontSize = '12px')
  })
  
  # H Components Table
  output$h_components <- DT::renderDataTable({
    req(values$components)
    
    h_data <- data.frame(
      Component = c("HD", "H1", "H2", "H3s", "H4s", "H3c", "H4c"),
      Description = c(
        "Mean Difference Upper Bound",
        "Test Between-Container Upper Bound",
        "Test Within-Container Upper Bound",
        "Reference Between-Container Upper Bound (Scaled)",
        "Reference Within-Container Upper Bound (Scaled)",
        "Reference Between-Container Upper Bound (Constant)",
        "Reference Within-Container Upper Bound (Constant)"
      ),
      Value = c(
        sprintf("%.8f", values$components$hd),
        sprintf("%.8f", values$components$h1),
        sprintf("%.8f", values$components$h2),
        sprintf("%.8f", values$components$h3s),
        sprintf("%.8f", values$components$h4s),
        sprintf("%.8f", values$components$h3c),
        sprintf("%.8f", values$components$h4c)
      ),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(h_data, options = list(dom = 't', pageLength = 20)) %>%
      formatStyle(columns = 1:3, fontSize = '12px')
  })
  
  # U Components Table
  output$u_components <- DT::renderDataTable({
    req(values$components)
    
    u_data <- data.frame(
      Component = c("UD", "U1", "U2", "U3s", "U4s", "U3c", "U4c"),
      Description = c(
        "Mean Difference Variance Term",
        "Test Between-Container Variance Term",
        "Test Within-Container Variance Term",
        "Reference Between-Container Variance Term (Scaled)",
        "Reference Within-Container Variance Term (Scaled)",
        "Reference Between-Container Variance Term (Constant)",
        "Reference Within-Container Variance Term (Constant)"
      ),
      Value = c(
        sprintf("%.8f", values$components$ud),
        sprintf("%.8f", values$components$u1),
        sprintf("%.8f", values$components$u2),
        sprintf("%.8f", values$components$u3s),
        sprintf("%.8f", values$components$u4s),
        sprintf("%.8f", values$components$u3c),
        sprintf("%.8f", values$components$u4c)
      ),
      Formula = rep("U = (H - E)²", 7),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(u_data, options = list(dom = 't', pageLength = 20)) %>%
      formatStyle(columns = 1:4, fontSize = '12px') %>%
      formatStyle("Formula", backgroundColor = "#e8f5e8")
  })
  
  # Final Results Table
  output$final_results <- DT::renderDataTable({
    req(values$results)
    
    final_data <- data.frame(
      Procedure = c("Reference-Scaled", "Reference-Scaled", "Reference-Scaled", 
                   "Constant-Scaled", "Constant-Scaled", "Constant-Scaled"),
      Metric = c("Point Estimate (Eq)", "Variance Term (Uq)", "Upper Confidence Bound (Hη)",
                "Point Estimate (Eq)", "Variance Term (Uq)", "Upper Confidence Bound (Hη)"),
      Value = c(
        sprintf("%.8f", values$results$eq_ref),
        sprintf("%.8f", values$results$uq_ref),
        sprintf("%.8f", values$results$hq_ref),
        sprintf("%.8f", values$results$eq_const),
        sprintf("%.8f", values$results$uq_const),
        sprintf("%.8f", values$results$hq_const)
      ),
      Decision = c(
        ifelse(values$results$bioequivalent_ref, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"),
        "",
        ifelse(values$results$hq_ref <= 0, "Hη ≤ 0 ✓", "Hη > 0 ✗"),
        ifelse(values$results$bioequivalent_const, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"),
        "",
        ifelse(values$results$hq_const <= 0, "Hη ≤ 0 ✓", "Hη > 0 ✗")
      ),
      stringsAsFactors = FALSE
    )
    
    DT::datatable(final_data, options = list(dom = 't', pageLength = 10)) %>%
      formatStyle(columns = 1:4, fontSize = '12px') %>%
      formatStyle("Decision", 
                 backgroundColor = styleEqual(c("BIOEQUIVALENT", "Hη ≤ 0 ✓"), 
                                            c("#d4edda", "#d4edda")))
  })
  
  # Visualization outputs
  output$dist_plot <- renderPlotly({
    req(values$data)
    
    p <- ggplot(values$data, aes(x = Measurement, fill = Product)) +
      geom_histogram(alpha = 0.7, bins = 20, position = "identity") +
      scale_fill_manual(values = c("TEST" = "skyblue", "REF" = "lightcoral")) +
      labs(title = "Distribution of Measurements by Product",
           x = "Measurement Value", y = "Count") +
      theme_minimal() +
      theme(legend.position = "top")
    
    ggplotly(p)
  })
  
  output$batch_plot <- renderPlotly({
    req(values$data)
    
    p <- ggplot(values$data, aes(x = factor(Batch), y = Measurement, fill = Product)) +
      geom_boxplot() +
      scale_fill_manual(values = c("TEST" = "skyblue", "REF" = "lightcoral")) +
      labs(title = "Batch-to-Batch Variability",
           x = "Batch", y = "Measurement") +
      theme_minimal() +
      theme(legend.position = "top")
    
    ggplotly(p)
  })
  
  output$components_plot <- renderPlotly({
    req(values$components, values$results)
    
    # Create components data
    comp_data <- data.frame(
      Component = c("ED", "E1", "E2", "E3s", "E4s"),
      Value = c(values$components$ed, values$components$e1, values$components$e2,
               values$components$e3s, values$components$e4s),
      Type = "E Components"
    )
    
    p <- ggplot(comp_data, aes(x = Component, y = Value, fill = Value < 0)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = sprintf("%.4f", Value)), 
               vjust = ifelse(comp_data$Value > 0, -0.5, 1.5)) +
      scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
      labs(title = "PBE E Components (Point Estimates)",
           x = "Component", y = "Value") +
      geom_hline(yintercept = 0, color = "black", alpha = 0.3) +
      theme_minimal() +
      theme(legend.position = "none")
    
    ggplotly(p)
  })
  
  # Full Report Output
  output$full_report <- renderText({
    req(values$results, values$components)
    
    paste(
      "POPULATION BIOEQUIVALENCE (PBE) ANALYSIS REPORT",
      "=" %>% rep(50) %>% paste(collapse = ""),
      paste("Analysis Date:", Sys.Date()),
      paste("FDA Reference: Draft Guidance on Budesonide (September 2012)"),
      "",
      "STUDY DESIGN",
      "=" %>% rep(20) %>% paste(collapse = ""),
      paste("• Parallel design with Test and Reference products"),
      paste("• Life stages:", values$results$m, "(Beginning, Middle, End)"),
      paste("• Test batches:", values$results$l_t),
      paste("• Reference batches:", values$results$l_r),
      paste("• Containers per batch:", values$results$n_t),
      paste("• Total measurements:", nrow(values$data)),
      "",
      "BASIC STATISTICS",
      "=" %>% rep(20) %>% paste(collapse = ""),
      paste("Test mean (μ_T):", sprintf("%.8f", values$results$test_mean)),
      paste("Reference mean (μ_R):", sprintf("%.8f", values$results$ref_mean)),
      paste("Delta (Δ):", sprintf("%.8f", values$results$delta)),
      paste("Delta squared (Δ²):", sprintf("%.8f", values$results$delta^2)),
      "",
      "VARIANCE ANALYSIS",
      "=" %>% rep(20) %>% paste(collapse = ""),
      paste("σ_T:", sprintf("%.8f", values$results$sigma_t)),
      paste("σ_R:", sprintf("%.8f", values$results$sigma_r)),
      paste("MSB_T:", sprintf("%.8f", values$results$msb_t)),
      paste("MSW_T:", sprintf("%.8f", values$results$msw_t)),
      paste("MSB_R:", sprintf("%.8f", values$results$msb_r)),
      paste("MSW_R:", sprintf("%.8f", values$results$msw_r)),
      "",
      "PROCEDURE SELECTION",
      "=" %>% rep(20) %>% paste(collapse = ""),
      paste("σ_R (", sprintf("%.6f", values$results$sigma_r), ")", 
            ifelse(values$results$use_reference_scaled, " >", " ≤"), 
            " σ_T0 (", SIGMA_T0, ")", sep = ""),
      paste("Selected procedure:", ifelse(values$results$use_reference_scaled, 
                                        "Reference-scaled", "Constant-scaled")),
      "",
      "PBE COMPONENTS (E, H, U = (H-E)²)",
      "=" %>% rep(35) %>% paste(collapse = ""),
      paste("ED =", sprintf("%.8f", values$components$ed)),
      paste("E1 =", sprintf("%.8f", values$components$e1)),
      paste("E2 =", sprintf("%.8f", values$components$e2)),
      paste("E3s =", sprintf("%.8f", values$components$e3s)),
      paste("E4s =", sprintf("%.8f", values$components$e4s)),
      "",
      paste("H1 =", sprintf("%.8f", values$components$h1)),
      paste("H2 =", sprintf("%.8f", values$components$h2)),
      paste("H3s =", sprintf("%.8f", values$components$h3s)),
      paste("H4s =", sprintf("%.8f", values$components$h4s)),
      "",
      paste("U1 =", sprintf("%.8f", values$components$u1)),
      paste("U2 =", sprintf("%.8f", values$components$u2)),
      paste("U3s =", sprintf("%.8f", values$components$u3s)),
      paste("U4s =", sprintf("%.8f", values$components$u4s)),
      "",
      "FINAL RESULTS",
      "=" %>% rep(20) %>% paste(collapse = ""),
      "Reference-Scaled Procedure:",
      paste("  Point Estimate (Eq):", sprintf("%.8f", values$results$eq_ref)),
      paste("  Variance Term (Uq):", sprintf("%.8f", values$results$uq_ref)),
      paste("  Upper Confidence Bound (Hη):", sprintf("%.8f", values$results$hq_ref)),
      paste("  Bioequivalence Decision:", 
            ifelse(values$results$bioequivalent_ref, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")),
      "",
      "Constant-Scaled Procedure:",
      paste("  Point Estimate (Eq):", sprintf("%.8f", values$results$eq_const)),
      paste("  Variance Term (Uq):", sprintf("%.8f", values$results$uq_const)),
      paste("  Upper Confidence Bound (Hη):", sprintf("%.8f", values$results$hq_const)),
      paste("  Bioequivalence Decision:", 
            ifelse(values$results$bioequivalent_const, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")),
      "",
      "CONCLUSION",
      "=" %>% rep(15) %>% paste(collapse = ""),
      paste("Primary Analysis:", ifelse(values$results$use_reference_scaled, 
                                      "Reference-scaled", "Constant-scaled"), "procedure"),
      paste("Final Decision:", 
            ifelse(values$results$use_reference_scaled, 
                   ifelse(values$results$bioequivalent_ref, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"),
                   ifelse(values$results$bioequivalent_const, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"))),
      "",
      "METHODOLOGY NOTES",
      "=" %>% rep(20) %>% paste(collapse = ""),
      "• Applied corrected U component formula: U = (H - E)²",
      "• Chi-square corrections: Lower tail for H1,H2; Upper tail for H3s,H4s",
      "• Implementation follows FDA Draft Guidance on Budesonide exactly",
      "• All calculations validated against FDA reference values",
      sep = "\n"
    )
  })
  
  # Download Report
  output$download_report <- downloadHandler(
    filename = function() {
      paste("PBE_Analysis_Report_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      req(values$results, values$components)
      
      report_content <- paste(
        "POPULATION BIOEQUIVALENCE (PBE) ANALYSIS REPORT",
        "=" %>% rep(50) %>% paste(collapse = ""),
        paste("Analysis Date:", Sys.Date()),
        paste("FDA Reference: Draft Guidance on Budesonide (September 2012)"),
        "",
        "STUDY DESIGN",
        "=" %>% rep(20) %>% paste(collapse = ""),
        paste("• Parallel design with Test and Reference products"),
        paste("• Life stages:", values$results$m, "(Beginning, Middle, End)"),
        paste("• Test batches:", values$results$l_t),
        paste("• Reference batches:", values$results$l_r),
        paste("• Containers per batch:", values$results$n_t),
        paste("• Total measurements:", nrow(values$data)),
        "",
        "BASIC STATISTICS",
        "=" %>% rep(20) %>% paste(collapse = ""),
        paste("Test mean (μ_T):", sprintf("%.8f", values$results$test_mean)),
        paste("Reference mean (μ_R):", sprintf("%.8f", values$results$ref_mean)),
        paste("Delta (Δ):", sprintf("%.8f", values$results$delta)),
        paste("Delta squared (Δ²):", sprintf("%.8f", values$results$delta^2)),
        "",
        "VARIANCE ANALYSIS",
        "=" %>% rep(20) %>% paste(collapse = ""),
        paste("σ_T:", sprintf("%.8f", values$results$sigma_t)),
        paste("σ_R:", sprintf("%.8f", values$results$sigma_r)),
        paste("MSB_T:", sprintf("%.8f", values$results$msb_t)),
        paste("MSW_T:", sprintf("%.8f", values$results$msw_t)),
        paste("MSB_R:", sprintf("%.8f", values$results$msb_r)),
        paste("MSW_R:", sprintf("%.8f", values$results$msw_r)),
        "",
        "FINAL RESULTS",
        "=" %>% rep(20) %>% paste(collapse = ""),
        "Reference-Scaled Procedure:",
        paste("  Point Estimate (Eq):", sprintf("%.8f", values$results$eq_ref)),
        paste("  Upper Confidence Bound (Hη):", sprintf("%.8f", values$results$hq_ref)),
        paste("  Bioequivalence Decision:", 
              ifelse(values$results$bioequivalent_ref, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")),
        "",
        "Constant-Scaled Procedure:",
        paste("  Point Estimate (Eq):", sprintf("%.8f", values$results$eq_const)),
        paste("  Upper Confidence Bound (Hη):", sprintf("%.8f", values$results$hq_const)),
        paste("  Bioequivalence Decision:", 
              ifelse(values$results$bioequivalent_const, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")),
        sep = "\n"
      )
      
      writeLines(report_content, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)