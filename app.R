# Population Bioequivalence (PBE) Analysis Shiny App
# Based on FDA Draft Guidance on Budesonide (September 2012)

library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(readr)
library(readxl)
library(openxlsx)
library(haven)
library(dplyr)
library(knitr)
library(kableExtra)
library(rmarkdown)
library(stringr)

# Define regulatory constants
SIGMA_T0 <- 0.1
THETA_P <- 2.0891

# Function to automatically detect CSV separator
detect_separator <- function(file_path, n_lines = 5) {
  # Read first few lines to detect separator
  sample_lines <- readLines(file_path, n = n_lines)
  
  # Count occurrences of common separators
  sep_counts <- list(
    "," = sum(str_count(sample_lines, ",")),
    ";" = sum(str_count(sample_lines, ";")),
    "\t" = sum(str_count(sample_lines, "\t"))
  )
  
  # Return the separator with the highest count
  best_sep <- names(sep_counts)[which.max(sep_counts)]
  
  # If no clear winner, default to comma
  if (sep_counts[[best_sep]] == 0) {
    return(",")
  }
  
  return(best_sep)
}

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
  
  # HD Component (t-distribution) - Fixed formula based on FDA PSG page 8
  # The FDA formula shows: HD = (|Δ| + t_{1-α,df} × √(variance))^2
  # where variance = MSB_T/(n_T × l_T × m) + MSB_R/(n_R × l_R × m)
  # and df = (n_T × l_T - 1) + (n_R × l_R - 1)
  df_pooled <- (l_t * n_t - 1) + (l_r * n_r - 1)
  t_critical <- qt(1 - alpha, df_pooled)  # One-tailed test: 1-α = 0.95
  # For budesonide example: df = (10×3-1) + (10×3-1) = 58, α = 0.05, so t_{0.95,58}
  hd_variance <- msb_t/(n_t * l_t * m) + msb_r/(n_r * l_r * m)
  hd <- (abs(delta) + t_critical * sqrt(hd_variance))^2
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
  ref_components$u3 <- u3s
  
  # E4s Component (Reference MSW - scaled)
  e4s <- -(1 + theta_p) * (m - 1) * msw_r / m
  ref_components$e4s <- e4s
  
  # H4s Component (Chi-square) - Upper tail for reference product
  df4 <- l_r * n_r * (m - 1)
  chi2_critical_4 <- qchisq(1 - alpha, df4)  # Upper tail for H3s, H4s
  h4s <- l_r * n_r * (m - 1) * e4s / chi2_critical_4
  u4s <- (h4s - e4s)^2  # CORRECTED: U = (H-E)²
  
  ref_components$h4s <- h4s
  ref_components$u4 <- u4s
  
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
ui <- page_sidebar(
  title = "Population Bioequivalence (PBE) Analysis",
  theme = bs_theme(version = 5, bootswatch = "minty"),
  
  # Sidebar with Upload Data functionality
  sidebar = sidebar(
    width = 350,
    
    # Data Upload Section
    h4("Data Upload"),
    fileInput("file", "Choose CSV, Excel, or XPT File",
             accept = c(".csv", ".xlsx", ".xls", ".xpt")),
    
    tags$hr(),
    
    # Column Mapping Section
    conditionalPanel(
      condition = "output.show_mapping",
      h4("Column Mapping"),
      p(tags$small("Map your columns to required parameters:")),
      
      div(style = "display: flex; align-items: center; margin-bottom: 5px;",
          "Batch:",
          tags$span(title = "Batch numbers or identifiers (e.g., 1, 2, 3 or B1, B2, B3)",
                   style = "margin-left: 8px; color: #007bff; cursor: default; user-select: none;",
                   HTML("&#9432;"))
      ),
      selectInput("map_batch", label = NULL, choices = NULL, width = "100%"),
      
      div(style = "display: flex; align-items: center; margin-bottom: 5px;",
          "Container:",
          tags$span(title = "Container/bottle IDs within each batch (e.g., 1-30, C1-C30)",
                   style = "margin-left: 8px; color: #007bff; cursor: default; user-select: none;",
                   HTML("&#9432;"))
      ),
      selectInput("map_container", label = NULL, choices = NULL, width = "100%"),
      
      div(style = "display: flex; align-items: center; margin-bottom: 5px;",
          "Stage:",
          tags$span(title = "Life stages - accepts B/Beginning, M/Middle, E/End (at least one required)",
                   style = "margin-left: 8px; color: #007bff; cursor: default; user-select: none;",
                   HTML("&#9432;"))
      ),
      selectInput("map_stage", label = NULL, choices = NULL, width = "100%"),
      
      div(style = "display: flex; align-items: center; margin-bottom: 5px;",
          "Product:",
          tags$span(title = "Product type - accepts TEST/REF, T/R, Test/Reference",
                   style = "margin-left: 8px; color: #007bff; cursor: default; user-select: none;",
                   HTML("&#9432;"))
      ),
      selectInput("map_product", label = NULL, choices = NULL, width = "100%"),
      
      div(style = "display: flex; align-items: center; margin-bottom: 5px;",
          "Measurement:",
          tags$span(title = "Numeric measurement values (any decimal format)",
                   style = "margin-left: 8px; color: #007bff; cursor: default; user-select: none;",
                   HTML("&#9432;"))
      ),
      selectInput("map_measurement", label = NULL, choices = NULL, width = "100%"),
      
      br()
    ),
    
    tags$hr(),
    
    # Analysis Button
    actionButton("analyze", "Run PBE Analysis", 
                class = "btn-success btn-lg", width = "100%"),
    
    br(), br()
  ),
  
  # Main Panel with Report
  card(
    card_header(
      div(
        style = "display: flex; justify-content: space-between; align-items: center;",
        h3("PBE Analysis Report - FDA PSG Format", style = "margin: 0;"),
        downloadButton("download_report", "Download Report", 
                      class = "btn-primary")
      )
    ),
    uiOutput("html_report")
  )
)

# Server
server <- function(input, output, session) {
  # Reactive values
  values <- reactiveValues(
    raw_data = NULL,      # Original uploaded data
    data = NULL,          # Mapped and processed data
    results = NULL,
    components = NULL,
    show_mapping = FALSE  # Control column mapping visibility
  )
  
  # File upload and column mapping setup
  observeEvent(input$file, {
    req(input$file)
    
    tryCatch({
      # Determine file type based on extension
      file_ext <- tools::file_ext(input$file$name)
      
      if (file_ext %in% c("xlsx", "xls")) {
        # Read Excel file
        values$raw_data <- read_excel(input$file$datapath, col_names = TRUE)
      } else if (file_ext == "xpt") {
        # Read SAS Transport file
        values$raw_data <- read_xpt(input$file$datapath)
      } else {
        # Read CSV file with automatic separator detection
        detected_sep <- detect_separator(input$file$datapath)
        values$raw_data <- read.csv(input$file$datapath,
                                   header = TRUE,
                                   sep = detected_sep)
      }
      
      # Get column names for mapping
      col_names <- names(values$raw_data)
      col_choices <- c("Select column..." = "", setNames(col_names, col_names))
      
      # Update column mapping dropdowns
      updateSelectInput(session, "map_batch", choices = col_choices, 
                       selected = ifelse("Batch" %in% col_names, "Batch", ""))
      updateSelectInput(session, "map_container", choices = col_choices,
                       selected = ifelse("Container" %in% col_names, "Container", ""))
      updateSelectInput(session, "map_stage", choices = col_choices,
                       selected = ifelse("Stage" %in% col_names, "Stage", ""))
      updateSelectInput(session, "map_product", choices = col_choices,
                       selected = ifelse("Product" %in% col_names, "Product", ""))
      updateSelectInput(session, "map_measurement", choices = col_choices,
                       selected = ifelse("Measurement" %in% col_names, "Measurement", ""))
      
      # Show column mapping panel
      values$show_mapping <- TRUE
      
      # Clear previous mapped data
      values$data <- NULL
      
      showNotification("File uploaded successfully! Please map the columns below.", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error", duration = 10)
      values$raw_data <- NULL
      values$show_mapping <- FALSE
    })
  })
  
  
  # Output to control column mapping visibility
  output$show_mapping <- reactive({
    values$show_mapping
  })
  outputOptions(output, "show_mapping", suspendWhenHidden = FALSE)
  
  
  # PBE Analysis
  observeEvent(input$analyze, {
    # First, handle column mapping if raw data exists but mapped data doesn't
    if (!is.null(values$raw_data) && is.null(values$data)) {
      # Check if all mappings are selected
      mappings <- list(
        Batch = input$map_batch,
        Container = input$map_container,
        Stage = input$map_stage,
        Product = input$map_product,
        Measurement = input$map_measurement
      )
      
      # Validate that all columns are mapped
      empty_mappings <- sapply(mappings, function(x) is.null(x) || x == "")
      if (any(empty_mappings)) {
        missing_params <- names(mappings)[empty_mappings]
        showNotification(paste("Please map all required columns. Missing:", paste(missing_params, collapse = ", ")), 
                        type = "error", duration = 10)
        return()
      }
      
      # Check for duplicate mappings
      selected_cols <- unlist(mappings)
      if (length(selected_cols) != length(unique(selected_cols))) {
        showNotification("Error: Each column can only be mapped once. Please select different columns for each parameter.", 
                        type = "error", duration = 10)
        return()
      }
      
      tryCatch({
        # Create mapped data frame with standardization
        raw_product <- values$raw_data[[input$map_product]]
        raw_stage <- values$raw_data[[input$map_stage]]
        
        # Standardize Product column
        standardized_product <- toupper(as.character(raw_product))
        standardized_product <- case_when(
          standardized_product %in% c("TEST", "T") ~ "TEST",
          standardized_product %in% c("REF", "REFERENCE", "R") ~ "REF",
          TRUE ~ standardized_product
        )
        
        # Standardize Stage column  
        standardized_stage <- toupper(as.character(raw_stage))
        standardized_stage <- case_when(
          standardized_stage %in% c("B", "BEGINNING", "BEGIN") ~ "B",
          standardized_stage %in% c("M", "MIDDLE", "MID") ~ "M", 
          standardized_stage %in% c("E", "END", "ENDING") ~ "E",
          TRUE ~ standardized_stage
        )
        
        values$data <- data.frame(
          Batch = values$raw_data[[input$map_batch]],
          Container = values$raw_data[[input$map_container]],
          Stage = standardized_stage,
          Product = standardized_product,
          Measurement = values$raw_data[[input$map_measurement]]
        )
        
        # Validate mapped data
        # Check for missing values
        if (any(is.na(values$data))) {
          showNotification("Warning: Data contains missing values. Please check your data.", 
                          type = "warning", duration = 10)
        }
        
        # Validate products
        products <- unique(values$data$Product)
        if (!all(c("TEST", "REF") %in% products)) {
          showNotification(paste("Error: Product column must contain both TEST and REF products. Found:", paste(products, collapse = ", ")), 
                          type = "error", duration = 10)
          values$data <- NULL
          return()
        }
        
        # Validate stages - just check that there's at least one stage
        stages <- unique(values$data$Stage)
        valid_stages <- c("B", "M", "E")
        if (!any(stages %in% valid_stages)) {
          showNotification(paste("Error: Stage column must contain at least one valid stage (B, M, or E). Found:", paste(stages, collapse = ", ")), 
                          type = "error", duration = 10)
          values$data <- NULL
          return()
        }
        
        # Validate measurement column is numeric
        if (!is.numeric(values$data$Measurement)) {
          # Try to convert to numeric
          values$data$Measurement <- as.numeric(values$data$Measurement)
          if (any(is.na(values$data$Measurement))) {
            showNotification("Error: Measurement column must contain numeric values", 
                            type = "error", duration = 10)
            values$data <- NULL
            return()
          }
        }
        
      }, error = function(e) {
        showNotification(paste("Error applying column mapping:", e$message), type = "error", duration = 10)
        values$data <- NULL
        return()
      })
    }
    
    # Now proceed with analysis
    req(values$data)
    
    tryCatch({
      # Separate test and reference data
      test_df <- filter(values$data, Product == 'TEST')
      ref_df <- filter(values$data, Product == 'REF')
      
      # Basic statistics - PBE analysis uses arithmetic means for delta calculation
      test_mean <- mean(test_df$Measurement)
      ref_mean <- mean(ref_df$Measurement)
      delta <- test_mean - ref_mean
      
      # Calculate geometric means properly
      test_log_mean <- mean(log(test_df$Measurement))
      ref_log_mean <- mean(log(ref_df$Measurement))
      test_geomean <- exp(test_log_mean)
      ref_geomean <- exp(ref_log_mean)
      
      # Study design parameters
      m <- length(unique(values$data$Stage))
      l_t <- length(unique(test_df$Batch))
      l_r <- length(unique(ref_df$Batch))
      n_t <- length(unique(test_df$Container)) / l_t
      n_r <- length(unique(ref_df$Container)) / l_r
      
      # Calculate MSB and MSW on original (raw) measurement data
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
      
      # Automatically determine PSG method based on μ_T and μ_R
      # Calculate μ_T and μ_R (arithmetic means)
      mu_t <- test_mean
      mu_r <- ref_mean
      
      # Determine PSG method based on PSG guidance criteria:
      # Traditional PBE (includes ED/UD): when μ_T > μ_R
      # Modified one-sided PBE (excludes ED/UD): when μ_T < μ_R
      mean_ratio <- mu_t / mu_r
      psg_method <- ifelse(mu_t < mu_r, "fluticasone", "budesonide")
      
      # Handle different PSG methods
      if (psg_method == "fluticasone") {
        # Fluticasone Propionate - EXCLUDE ED and UD terms
        # Reference-scaled procedure (without ED/UD)
        eq_ref <- all_components$e1 + all_components$e2 + 
                  all_components$e3s + all_components$e4s
        uq_ref <- all_components$u1 + all_components$u2 + 
                  all_components$u3 + all_components$u4
        hq_ref <- eq_ref + sqrt(max(0, uq_ref))
        bioequivalent_ref <- hq_ref <= 0
        
        # Constant-scaled procedure (without ED/UD)
        eq_const <- all_components$e1 + all_components$e2 + 
                    all_components$e3c + all_components$e4c - THETA_P * SIGMA_T0^2
        uq_const <- all_components$u1 + all_components$u2 + 
                    all_components$u3c + all_components$u4c
        hq_const <- eq_const + sqrt(max(0, uq_const))
        bioequivalent_const <- hq_const <= 0
        
        fluticasone_results <- list(
          procedure_used = "Fluticasone Propionate PSG - Excludes ED/UD terms",
          excludes_ed_ud = TRUE
        )
        
      } else {
        # Standard Budesonide PBE - INCLUDE ED and UD terms
        # Reference-scaled procedure (with ED/UD)
        eq_ref <- all_components$ed + all_components$e1 + all_components$e2 + 
                  all_components$e3s + all_components$e4s
        uq_ref <- all_components$ud + all_components$u1 + all_components$u2 + 
                  all_components$u3 + all_components$u4
        hq_ref <- eq_ref + sqrt(max(0, uq_ref))
        bioequivalent_ref <- hq_ref <= 0
        
        # Constant-scaled procedure (with ED/UD)
        eq_const <- all_components$ed + all_components$e1 + all_components$e2 + 
                    all_components$e3c + all_components$e4c - THETA_P * SIGMA_T0^2
        uq_const <- all_components$ud + all_components$u1 + all_components$u2 + 
                    all_components$u3c + all_components$u4c
        hq_const <- eq_const + sqrt(max(0, uq_const))
        bioequivalent_const <- hq_const <= 0
        
        fluticasone_results <- NULL
      }
      
      # Store results
      values$results <- list(
        test_mean = test_mean,
        ref_mean = ref_mean,
        mu_t = mu_t,
        mu_r = mu_r,
        mean_ratio = mean_ratio,
        test_log_mean = test_log_mean,
        ref_log_mean = ref_log_mean,
        test_geomean = test_geomean,
        ref_geomean = ref_geomean,
        delta = delta,
        sigma_t = sigma_t,
        sigma_r = sigma_r,
        msb_t = msb_t, msw_t = msw_t,
        msb_r = msb_r, msw_r = msw_r,
        m = m, l_t = l_t, l_r = l_r, n_t = n_t, n_r = n_r,
        use_reference_scaled = use_reference_scaled,
        eq_ref = eq_ref, uq_ref = uq_ref, hq_ref = hq_ref, bioequivalent_ref = bioequivalent_ref,
        eq_const = eq_const, uq_const = uq_const, hq_const = hq_const, bioequivalent_const = bioequivalent_const,
        psg_method = psg_method,
        fluticasone_results = fluticasone_results
      )
      
      values$components <- all_components
      
      showNotification("PBE Analysis completed successfully!", type = "message")
      # Switch to results tab after analysis
      # updateTabItems(session, "tabs", "results")
      
    }, error = function(e) {
      showNotification(paste("Error in analysis:", e$message), type = "error", duration = 10)
    })
  })
  
  # HTML Report Generation
  output$html_report <- renderUI({
    req(values$results, values$components, values$data)
    
    tryCatch({
      # Generate HTML directly in R
      results <- values$results
      components <- values$components
      data_summary <- values$data
      
      # Create HTML content
      div(
        # Header
        div(
          style = "text-align: center; margin-bottom: 30px;",
          h2("Population Bioequivalence (PBE) Analysis Report"),
          h4("FDA Draft Guidance Format"),
          p(paste("Generated:", Sys.Date()))
        ),
        
        # Study Summary
        h3("Study Summary"),
        p(strong("PSG Method:"), ifelse(results$psg_method == "fluticasone", "Modified PBE (Fluticasone-type)", "Standard PBE (Budesonide-type)"), 
          " - ", em("(automatically determined)")),
        p(strong("Method Selection Criteria:"), paste("μ_T =", sprintf("%.4f", results$mu_t), ", μ_R =", sprintf("%.4f", results$mu_r), 
                                                     ifelse(results$psg_method == "fluticasone", "(μ_T < μ_R, excludes ED/UD)", "(μ_T ≥ μ_R, includes ED/UD)"))),
        p(strong("Selected Procedure:"), ifelse(results$use_reference_scaled, "Reference-scaled", "Constant-scaled")),
        p(strong("Procedure Selection:"), paste("σ_R (", sprintf("%.6f", results$sigma_r), ")", ifelse(results$use_reference_scaled, " > ", " ≤ "), "σ_T0 (0.1)")),
        
        hr(),
        
        # Reference-Scaled Analysis Results
        h3("Reference-Scaled Analysis Results"),
        if(results$psg_method == "fluticasone") {
          # Fluticasone PSG - exclude ED/UD terms
          ref_data <- data.frame(
            col1 = c(
              sprintf("E1 = %.9f", components$e1), 
              sprintf("E2 = %.8f", components$e2),
              sprintf("E3s = %.9f", components$e3s),
              sprintf("E4s = %.9f", components$e4s),
              sprintf("Eq = %.8f", results$eq_ref)
            ),
            col2 = c(
              sprintf("H1 = %.9f", components$h1),
              sprintf("H2 = %.8f", components$h2), 
              sprintf("H3s = %.9f", components$h3s),
              sprintf("H4s = %.9f", components$h4s),
              ""
            ),
            col3 = c(
              sprintf("U1 = %.8f", components$u1),
              sprintf("U2 = %.9f", components$u2),
              sprintf("U3 = %.8f", components$u3),
              sprintf("U4 = %.8f", components$u4),
              sprintf("Uq = %.9f", results$uq_ref)
            ),
            col4 = c(
              "", "", "", "",
              sprintf("Hη = %.9f", results$hq_ref)
            )
          )
          
          ref_data %>%
            kable(format = "html", escape = FALSE, 
                  col.names = c("Point Estimates (Eq)", "Upper Bounds (Hq)", "Variance Terms (Uq)", "Final Result"),
                  caption = "Reference-Scaled Analysis") %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                          full_width = TRUE, font_size = 14) %>%
            row_spec(5, bold = TRUE, background = "#f0f8ff") %>%
            HTML()
            
        } else {
          # Budesonide PSG - include ED/UD terms
          ref_data <- data.frame(
            col1 = c(
              sprintf("ED = %.9f", components$ed),
              sprintf("E1 = %.9f", components$e1), 
              sprintf("E2 = %.8f", components$e2),
              sprintf("E3s = %.9f", components$e3s),
              sprintf("E4s = %.9f", components$e4s),
              sprintf("Eq = %.8f", results$eq_ref)
            ),
            col2 = c(
              sprintf("HD = %.9f", components$hd),
              sprintf("H1 = %.9f", components$h1),
              sprintf("H2 = %.8f", components$h2), 
              sprintf("H3s = %.9f", components$h3s),
              sprintf("H4s = %.9f", components$h4s),
              ""
            ),
            col3 = c(
              sprintf("UD = %.9f", components$ud),
              sprintf("U1 = %.8f", components$u1),
              sprintf("U2 = %.9f", components$u2),
              sprintf("U3 = %.8f", components$u3),
              sprintf("U4 = %.8f", components$u4),
              sprintf("Uq = %.9f", results$uq_ref)
            ),
            col4 = c(
              "", "", "", "", "",
              sprintf("Hη = %.9f", results$hq_ref)
            )
          )
          
          ref_data %>%
            kable(format = "html", escape = FALSE, 
                  col.names = c("Point Estimates (Eq)", "Upper Bounds (Hq)", "Variance Terms (Uq)", "Final Result"),
                  caption = "Reference-Scaled Analysis") %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                          full_width = TRUE, font_size = 14) %>%
            row_spec(6, bold = TRUE, background = "#f0f8ff") %>%
            HTML()
        },
        
        hr(),
        
        # Constant-Scaled Analysis Results
        h3("Constant-Scaled Analysis Results"),
        if(results$psg_method == "fluticasone") {
          # Fluticasone PSG - exclude ED/UD terms
          const_data <- data.frame(
            col1 = c(
              sprintf("E1 = %.9f", components$e1),
              sprintf("E2 = %.8f", components$e2),
              sprintf("E3c = %.9f", components$e3c),
              sprintf("E4c = %.8f", components$e4c),
              sprintf("Eq = %.9f", results$eq_const)
            ),
            col2 = c(
              sprintf("H1 = %.9f", components$h1),
              sprintf("H2 = %.8f", components$h2),
              sprintf("H3c = %.9f", components$h3c),
              sprintf("H4c = %.8f", components$h4c),
              ""
            ),
            col3 = c(
              sprintf("U1 = %.8f", components$u1),
              sprintf("U2 = %.9f", components$u2),
              sprintf("U3c = %.9f", components$u3c),
              sprintf("U4c = %.9f", components$u4c),
              sprintf("Uq = %.9f", results$uq_const)
            ),
            col4 = c(
              "", "", "", "",
              sprintf("Hη = %.9f", results$hq_const)
            )
          )
          
          const_data %>%
            kable(format = "html", escape = FALSE,
                  col.names = c("Point Estimates (Eq)", "Upper Bounds (Hq)", "Variance Terms (Uq)", "Final Result"),
                  caption = "Constant-Scaled Analysis") %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                          full_width = TRUE, font_size = 14) %>%
            row_spec(5, bold = TRUE, background = "#fff8f0") %>%
            HTML()
            
        } else {
          # Budesonide PSG - include ED/UD terms
          const_data <- data.frame(
            col1 = c(
              sprintf("ED = %.9f", components$ed),
              sprintf("E1 = %.9f", components$e1),
              sprintf("E2 = %.8f", components$e2),
              sprintf("E3c = %.9f", components$e3c),
              sprintf("E4c = %.8f", components$e4c),
              sprintf("Eq = %.9f", results$eq_const)
            ),
            col2 = c(
              sprintf("HD = %.9f", components$hd),
              sprintf("H1 = %.9f", components$h1),
              sprintf("H2 = %.8f", components$h2),
              sprintf("H3c = %.9f", components$h3c),
              sprintf("H4c = %.8f", components$h4c),
              ""
            ),
            col3 = c(
              sprintf("UD = %.9f", components$ud),
              sprintf("U1 = %.8f", components$u1),
              sprintf("U2 = %.9f", components$u2),
              sprintf("U3c = %.9f", components$u3c),
              sprintf("U4c = %.9f", components$u4c),
              sprintf("Uq = %.9f", results$uq_const)
            ),
            col4 = c(
              "", "", "", "", "",
              sprintf("Hη = %.9f", results$hq_const)  
            )
          )
          
          const_data %>%
            kable(format = "html", escape = FALSE,
                  col.names = c("Point Estimates (Eq)", "Upper Bounds (Hq)", "Variance Terms (Uq)", "Final Result"),
                  caption = "Constant-Scaled Analysis") %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                          full_width = TRUE, font_size = 14) %>%
            row_spec(6, bold = TRUE, background = "#fff8f0") %>%
            HTML()
        },
        
        hr(),
        
        # Table 2: Basic Summary Statistics  
        h3("Table 2: Basic Summary Statistics"),
        {
          basic_stats_table <- data.frame(
            Variable = c(
              "In Vitro Measurement"
            ),
            `Geometric Mean Test` = c(
              sprintf("%.6f", results$test_geomean)
            ),
            `Geometric Mean Reference` = c(
              sprintf("%.6f", results$ref_geomean)
            ),
            `Geometric Mean Ratio` = c(
              sprintf("%.6f", results$test_geomean / results$ref_geomean)
            ),
            `Standard Deviation SigmaT` = c(
              sprintf("%.6f", results$sigma_t)
            ),
            `Standard Deviation SigmaR` = c(
              sprintf("%.6f", results$sigma_r)
            ),
            `SigmaT/SigmaR Ratio` = c(
              sprintf("%.6f", results$sigma_t / results$sigma_r)
            ),
            check.names = FALSE
          )
          
          basic_stats_table %>%
            kable(format = "html", escape = FALSE,
                  caption = "Table 2. Basic Summary Statistics") %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                          full_width = TRUE, font_size = 14) %>%
            row_spec(1, background = "#f8f9fa") %>%
            HTML()
        },
        
        
        hr(),
        
        # Final Bioequivalence Results
        h3("Final Bioequivalence Results"),
        {
          final_data <- data.frame(
            Procedure = c("Reference-Scaled", "Constant-Scaled"),
            `Point Estimate (Eq)` = c(
              sprintf("%.9f", results$eq_ref),
              sprintf("%.9f", results$eq_const)
            ),
            `Variance Term (Uq)` = c(
              sprintf("%.9f", results$uq_ref),
              sprintf("%.9f", results$uq_const)
            ),
            `Upper Confidence Bound (Hη)` = c(
              sprintf("%.9f", results$hq_ref),
              sprintf("%.9f", results$hq_const)
            ),
            Decision = c(
              ifelse(results$bioequivalent_ref, "BIOEQUIVALENT", "NOT BIOEQUIVALENT"),
              ifelse(results$bioequivalent_const, "BIOEQUIVALENT", "NOT BIOEQUIVALENT")
            ),
            check.names = FALSE
          )
          
          final_data %>%
            kable(format = "html", escape = FALSE,
                  caption = "Final Bioequivalence Results") %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                          full_width = TRUE, font_size = 14) %>%
            column_spec(5, bold = TRUE) %>%
            row_spec(which(results$use_reference_scaled), background = "#d4edda") %>%
            HTML()
        },
        
        hr(),
        
        # Appendix
        h3("Appendix"),
        h4("Table 1: Complete Individual Data of In Vitro Tests"),
        {
          full_data <- data_summary[, c("Batch", "Container", "Stage", "Product", "Measurement")]
          names(full_data) <- c("Batches", "Container", "Stage", "Product", "In vitro measurement (original data)")
          full_data$`In vitro measurement (original data)` <- sprintf("%.6f", full_data$`In vitro measurement (original data)`)
          
          full_data %>%
            kable(format = "html",
                  caption = "Table 1. Complete Individual Data of In Vitro Tests") %>%
            kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                          full_width = TRUE, font_size = 12) %>%
            scroll_box(width = "100%", height = "400px") %>%
            HTML()
        },
        
        # Data Visualizations
        h4("Data Visualizations"),
        h5("Data Distribution by Product"),
        renderPlot({
          ggplot(data_summary, aes(x = Measurement, fill = Product)) +
            geom_histogram(alpha = 0.7, bins = 20, position = "identity") +
            scale_fill_manual(values = c("TEST" = "skyblue", "REF" = "lightcoral")) +
            labs(title = "Distribution of Measurements by Product",
                 x = "Measurement Value", y = "Count") +
            theme_minimal() +
            theme(legend.position = "top",
                  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        }),
        
        h5("Batch-to-Batch Variability"),  
        renderPlot({
          ggplot(data_summary, aes(x = factor(Batch), y = Measurement, fill = Product)) +
            geom_boxplot() +
            scale_fill_manual(values = c("TEST" = "skyblue", "REF" = "lightcoral")) +
            labs(title = "Batch-to-Batch Variability",
                 x = "Batch", y = "Measurement") +
            theme_minimal() +
            theme(legend.position = "top",
                  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        }),
        
        hr(),
        p(strong("Report Generated:"), Sys.time())
      )
      
    }, error = function(e) {
      # Fallback in case of error
      div(
        h3("Report Generation Error"),
        p("There was an error generating the report. Please ensure all calculations have completed."),
        p(paste("Error:", e$message))
      )
    })
  })
  
  
  # Download Report as HTML
  output$download_report <- downloadHandler(
    filename = function() {
      paste("PBE_Analysis_Report_", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      req(values$results, values$components, values$data)
      
      tryCatch({
        # Generate the same HTML content as displayed in the app
        results <- values$results
        components <- values$components
        data_summary <- values$data
        
        # Create complete HTML document
        html_content <- paste0(
          '<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Population Bioequivalence (PBE) Analysis Report</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
  <style>
    body { font-family: Arial, sans-serif; margin: 20px; }
    .table { margin: 20px 0; }
    .text-center { text-align: center; }
    .mb-4 { margin-bottom: 1.5rem; }
    .mt-4 { margin-top: 1.5rem; }
    hr { margin: 2rem 0; }
  </style>
</head>
<body>
  <div class="container">
    <div class="text-center mb-4">
      <h2>Population Bioequivalence (PBE) Analysis Report</h2>
      <h4>FDA Draft Guidance Format</h4>
      <p>Generated: ', Sys.Date(), '</p>
    </div>
    
    <h3>Study Summary</h3>
    <p><strong>PSG Method:</strong> ', ifelse(results$psg_method == "fluticasone", "Modified PBE (Fluticasone-type)", "Standard PBE (Budesonide-type)"), ' - <em>(automatically determined)</em></p>
    <p><strong>Method Selection Criteria:</strong> μ_T = ', sprintf("%.4f", results$mu_t), ', μ_R = ', sprintf("%.4f", results$mu_r), ifelse(results$psg_method == "fluticasone", " (μ_T < μ_R, excludes ED/UD)", " (μ_T ≥ μ_R, includes ED/UD)"), '</p>
    <p><strong>Selected Procedure:</strong> ', ifelse(results$use_reference_scaled, "Reference-scaled", "Constant-scaled"), '</p>
    <p><strong>Procedure Selection:</strong> σ_R (', sprintf("%.6f", results$sigma_r), ')', ifelse(results$use_reference_scaled, " > ", " ≤ "), 'σ_T0 (0.1)</p>
    
    <hr>
    
    <h3>Reference-Scaled Analysis Results</h3>',
    
    # Reference-scaled table
    if(results$psg_method == "fluticasone") {
      # Fluticasone PSG - exclude ED/UD terms
      ref_data <- data.frame(
        col1 = c(
          sprintf("E1 = %.9f", components$e1), 
          sprintf("E2 = %.8f", components$e2),
          sprintf("E3s = %.9f", components$e3s),
          sprintf("E4s = %.9f", components$e4s),
          sprintf("<strong>Eq = %.8f</strong>", results$eq_ref)
        ),
        col2 = c(
          sprintf("H1 = %.9f", components$h1),
          sprintf("H2 = %.8f", components$h2), 
          sprintf("H3s = %.9f", components$h3s),
          sprintf("H4s = %.9f", components$h4s),
          ""
        ),
        col3 = c(
          sprintf("U1 = %.8f", components$u1),
          sprintf("U2 = %.9f", components$u2),
          sprintf("U3 = %.8f", components$u3),
          sprintf("U4 = %.8f", components$u4),
          sprintf("<strong>Uq = %.9f</strong>", results$uq_ref)
        ),
        col4 = c(
          "", "", "", "",
          sprintf("<strong>Hη = %.9f</strong>", results$hq_ref)
        )
      )
      
      ref_table <- kable(ref_data, format = "html", escape = FALSE, 
                        col.names = c("Point Estimates (Eq)", "Upper Bounds (Hq)", "Variance Terms (Uq)", "Final Result"),
                        table.attr = 'class="table table-striped table-hover"') %>%
                   as.character()
                   
    } else {
      # Budesonide PSG - include ED/UD terms
      ref_data <- data.frame(
        col1 = c(
          sprintf("ED = %.9f", components$ed),
          sprintf("E1 = %.9f", components$e1), 
          sprintf("E2 = %.8f", components$e2),
          sprintf("E3s = %.9f", components$e3s),
          sprintf("E4s = %.9f", components$e4s),
          sprintf("<strong>Eq = %.8f</strong>", results$eq_ref)
        ),
        col2 = c(
          sprintf("HD = %.9f", components$hd),
          sprintf("H1 = %.9f", components$h1),
          sprintf("H2 = %.8f", components$h2), 
          sprintf("H3s = %.9f", components$h3s),
          sprintf("H4s = %.9f", components$h4s),
          ""
        ),
        col3 = c(
          sprintf("UD = %.9f", components$ud),
          sprintf("U1 = %.8f", components$u1),
          sprintf("U2 = %.9f", components$u2),
          sprintf("U3 = %.8f", components$u3),
          sprintf("U4 = %.8f", components$u4),
          sprintf("<strong>Uq = %.9f</strong>", results$uq_ref)
        ),
        col4 = c(
          "", "", "", "", "",
          sprintf("<strong>Hη = %.9f</strong>", results$hq_ref)
        )
      )
      
      ref_table <- kable(ref_data, format = "html", escape = FALSE, 
                        col.names = c("Point Estimates (Eq)", "Upper Bounds (Hq)", "Variance Terms (Uq)", "Final Result"),
                        table.attr = 'class="table table-striped table-hover"') %>%
                   as.character()
    },
    
    ref_table,
    
    '<hr>
    
    <h3>Constant-Scaled Analysis Results</h3>',
    
    # Constant-scaled table
    if(results$psg_method == "fluticasone") {
      # Fluticasone PSG - exclude ED/UD terms
      const_data <- data.frame(
        col1 = c(
          sprintf("E1 = %.9f", components$e1),
          sprintf("E2 = %.8f", components$e2),
          sprintf("E3c = %.9f", components$e3c),
          sprintf("E4c = %.8f", components$e4c),
          sprintf("<strong>Eq = %.9f</strong>", results$eq_const)
        ),
        col2 = c(
          sprintf("H1 = %.9f", components$h1),
          sprintf("H2 = %.8f", components$h2),
          sprintf("H3c = %.9f", components$h3c),
          sprintf("H4c = %.8f", components$h4c),
          ""
        ),
        col3 = c(
          sprintf("U1 = %.8f", components$u1),
          sprintf("U2 = %.9f", components$u2),
          sprintf("U3c = %.9f", components$u3c),
          sprintf("U4c = %.9f", components$u4c),
          sprintf("<strong>Uq = %.9f</strong>", results$uq_const)
        ),
        col4 = c(
          "", "", "", "",
          sprintf("<strong>Hη = %.9f</strong>", results$hq_const)
        )
      )
      
      const_table <- kable(const_data, format = "html", escape = FALSE,
                          col.names = c("Point Estimates (Eq)", "Upper Bounds (Hq)", "Variance Terms (Uq)", "Final Result"),
                          table.attr = 'class="table table-striped table-hover"') %>%
                     as.character()
                     
    } else {
      # Budesonide PSG - include ED/UD terms
      const_data <- data.frame(
        col1 = c(
          sprintf("ED = %.9f", components$ed),
          sprintf("E1 = %.9f", components$e1),
          sprintf("E2 = %.8f", components$e2),
          sprintf("E3c = %.9f", components$e3c),
          sprintf("E4c = %.8f", components$e4c),
          sprintf("<strong>Eq = %.9f</strong>", results$eq_const)
        ),
        col2 = c(
          sprintf("HD = %.9f", components$hd),
          sprintf("H1 = %.9f", components$h1),
          sprintf("H2 = %.8f", components$h2),
          sprintf("H3c = %.9f", components$h3c),
          sprintf("H4c = %.8f", components$h4c),
          ""
        ),
        col3 = c(
          sprintf("UD = %.9f", components$ud),
          sprintf("U1 = %.8f", components$u1),
          sprintf("U2 = %.9f", components$u2),
          sprintf("U3c = %.9f", components$u3c),
          sprintf("U4c = %.9f", components$u4c),
          sprintf("<strong>Uq = %.9f</strong>", results$uq_const)
        ),
        col4 = c(
          "", "", "", "", "",
          sprintf("<strong>Hη = %.9f</strong>", results$hq_const)
        )
      )
      
      const_table <- kable(const_data, format = "html", escape = FALSE,
                          col.names = c("Point Estimates (Eq)", "Upper Bounds (Hq)", "Variance Terms (Uq)", "Final Result"),
                          table.attr = 'class="table table-striped table-hover"') %>%
                     as.character()
    },
    
    const_table,
    
    '<hr>
    
    <h3>Table 2: Basic Summary Statistics</h3>',
    
    # Basic statistics table
    {
      basic_stats_table <- data.frame(
        Variable = c(
          "In Vitro Measurement"
        ),
        `Geometric Mean Test` = c(
          sprintf("%.6f", results$test_geomean)
        ),
        `Geometric Mean Reference` = c(
          sprintf("%.6f", results$ref_geomean)
        ),
        `Geometric Mean Ratio` = c(
          sprintf("%.6f", results$test_geomean / results$ref_geomean)
        ),
        `Standard Deviation SigmaT` = c(
          sprintf("%.6f", results$sigma_t)
        ),
        `Standard Deviation SigmaR` = c(
          sprintf("%.6f", results$sigma_r)
        ),
        `SigmaT/SigmaR Ratio` = c(
          sprintf("%.6f", results$sigma_t / results$sigma_r)
        ),
        check.names = FALSE
      )
      
      kable(basic_stats_table, format = "html", escape = FALSE,
            table.attr = 'class="table table-striped table-hover"') %>%
        as.character()
    },
    
    
    '<hr>
    
    <h3>Final Bioequivalence Results</h3>',
    
    # Final results table
    {
      final_data <- data.frame(
        Procedure = c("Reference-Scaled", "Constant-Scaled"),
        `Point Estimate (Eq)` = c(
          sprintf("%.9f", results$eq_ref),
          sprintf("%.9f", results$eq_const)
        ),
        `Variance Term (Uq)` = c(
          sprintf("%.9f", results$uq_ref),
          sprintf("%.9f", results$uq_const)
        ),
        `Upper Confidence Bound (Hη)` = c(
          sprintf("%.9f", results$hq_ref),
          sprintf("%.9f", results$hq_const)
        ),
        Decision = c(
          ifelse(results$bioequivalent_ref, "<strong>BIOEQUIVALENT</strong>", "<strong>NOT BIOEQUIVALENT</strong>"),
          ifelse(results$bioequivalent_const, "<strong>BIOEQUIVALENT</strong>", "<strong>NOT BIOEQUIVALENT</strong>")
        ),
        check.names = FALSE
      )
      
      kable(final_data, format = "html", escape = FALSE,
            table.attr = 'class="table table-striped table-hover"') %>%
        as.character()
    },
    
    '<hr>
    
    <h3>Appendix</h3>
    <h4>Table 1: Complete Individual Data of In Vitro Tests</h4>',
    
    # Individual data table
    {
      full_data <- data_summary[, c("Batch", "Container", "Stage", "Product", "Measurement")]
      names(full_data) <- c("Batches", "Container", "Stage", "Product", "In vitro measurement (original data)")
      full_data$`In vitro measurement (original data)` <- sprintf("%.6f", full_data$`In vitro measurement (original data)`)
      
      kable(full_data, format = "html",
            table.attr = 'class="table table-striped table-hover table-sm"') %>%
        as.character()
    },
    
    '<hr>
    <p><strong>Report Generated:</strong> ', as.character(Sys.time()), '</p>
  </div>
</body>
</html>'
        )
        
        # Write HTML to file
        writeLines(html_content, file, useBytes = TRUE)
        
      }, error = function(e) {
        # Fallback: create a simple HTML file if generation fails
        error_html <- paste0(
          '<!DOCTYPE html>
<html>
<head><title>PBE Analysis Report - Error</title></head>
<body>
  <h1>PBE Analysis Report</h1>
  <p>Generated: ', Sys.time(), '</p>
  <p>PSG Method: ', values$results$psg_method, '</p>
  <hr>
  <p><strong>Error generating report:</strong> ', e$message, '</p>
</body>
</html>'
        )
        writeLines(error_html, file, useBytes = TRUE)
      })
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
