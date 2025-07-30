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
  # The FDA formula shows: HD = (Δ + t_α/2,df × √(variance))^2
  # where variance = MSB_T/(n_T × l_T × m) + MSB_R/(n_R × l_R × m)
  # and df = (n_T × l_T - 1) + (n_R × l_R - 1)
  df_pooled <- (l_t * n_t - 1) + (l_r * n_r - 1)
  t_critical <- qt(1 - alpha, df_pooled)
  # For budesonide example: df = (10×3-1) + (10×3-1) = 58, α = 0.05
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
    
    # PSG Method Selection
    h4("PSG Method"),
    radioButtons("psg_method", 
                label = "Select Analysis Method:",
                choices = list(
                  "Budesonide PSG" = "budesonide",
                  "Fluticasone Propionate PSG" = "fluticasone"
                ),
                selected = "budesonide",
                width = "100%"),
    
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
          tags$span(title = "Life stages - accepts B/Beginning, M/Middle, E/End",
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
      
      br(),
      actionButton("apply_mapping", "Apply Mapping", 
                  class = "btn-primary btn-sm", width = "100%"),
      br(), br()
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
  
  # Handle column mapping application
  observeEvent(input$apply_mapping, {
    req(values$raw_data)
    
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
      
      # Validate stages
      stages <- unique(values$data$Stage)
      expected_stages <- c("B", "M", "E")
      if (!all(expected_stages %in% stages)) {
        showNotification(paste("Error: Stage column must contain B, M, E stages. Found:", paste(stages, collapse = ", ")), 
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
      
      showNotification("Column mapping applied successfully! Data is ready for analysis.", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error applying column mapping:", e$message), type = "error", duration = 10)
      values$data <- NULL
    })
  })
  
  # Output to control column mapping visibility
  output$show_mapping <- reactive({
    values$show_mapping
  })
  outputOptions(output, "show_mapping", suspendWhenHidden = FALSE)
  
  
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
      
      # Handle different PSG methods
      if (input$psg_method == "fluticasone") {
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
        delta = delta,
        sigma_t = sigma_t,
        sigma_r = sigma_r,
        msb_t = msb_t, msw_t = msw_t,
        msb_r = msb_r, msw_r = msw_r,
        m = m, l_t = l_t, l_r = l_r, n_t = n_t, n_r = n_r,
        use_reference_scaled = use_reference_scaled,
        eq_ref = eq_ref, uq_ref = uq_ref, hq_ref = hq_ref, bioequivalent_ref = bioequivalent_ref,
        eq_const = eq_const, uq_const = uq_const, hq_const = hq_const, bioequivalent_const = bioequivalent_const,
        psg_method = input$psg_method,
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
    
    # Create temporary file for rendered HTML
    temp_rmd <- tempfile(fileext = ".Rmd")
    temp_html <- tempfile(fileext = ".html")
    
    # Copy template to temp location
    template_path <- file.path("PBE_SHINY", "report_template.Rmd")
    if (!file.exists(template_path)) {
      template_path <- "report_template.Rmd"  # fallback
    }
    file.copy(template_path, temp_rmd)
    
    tryCatch({
      # Render the R Markdown with parameters
      rmarkdown::render(
        temp_rmd,
        output_file = temp_html,
        params = list(
          results = values$results,
          components = values$components,
          data_summary = values$data
        ),
        quiet = TRUE
      )
      
      # Read the rendered HTML
      html_content <- readLines(temp_html, warn = FALSE)
      html_string <- paste(html_content, collapse = "\n")
      
      # Clean up temp files
      unlink(c(temp_rmd, temp_html))
      
      # Return HTML content
      HTML(html_string)
      
    }, error = function(e) {
      # Fallback in case of rendering error
      div(
        h3("Report Generation Error"),
        p("There was an error generating the report. Please ensure all calculations have completed."),
        p(paste("Error:", e$message))
      )
    })
  })
  
  
  # Download Report as PDF
  output$download_report <- downloadHandler(
    filename = function() {
      paste("PBE_Analysis_Report_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      req(values$results, values$components, values$data)
      
      # Create temporary file for R Markdown
      temp_rmd <- tempfile(fileext = ".Rmd")
      temp_pdf <- tempfile(fileext = ".pdf")
      
      # Copy template to temp location
      template_path <- file.path("PBE_SHINY", "report_template.Rmd")
      if (!file.exists(template_path)) {
        template_path <- "report_template.Rmd"  # fallback
      }
      file.copy(template_path, temp_rmd)
      
      tryCatch({
        # Render the R Markdown to PDF with parameters
        rmarkdown::render(
          temp_rmd,
          output_format = "pdf_document",
          output_file = temp_pdf,
          params = list(
            results = values$results,
            components = values$components,
            data_summary = values$data
          ),
          quiet = TRUE
        )
        
        # Copy the generated PDF to the download file
        file.copy(temp_pdf, file)
        
        # Clean up temp files
        unlink(c(temp_rmd, temp_pdf))
        
      }, error = function(e) {
        # Fallback: create a simple text file if PDF generation fails
        report_content <- paste(
          "PBE Analysis Report",
          paste("Generated:", Sys.time()),
          paste("PSG Method:", values$results$psg_method),
          "",
          "PDF generation failed. Please check that you have LaTeX installed.",
          paste("Error:", e$message),
          sep = "\n"
        )
        writeLines(report_content, file)
      })
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
