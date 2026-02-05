library(shiny)
library(tidyverse)
library(plyr)
library(dplyr)
library(tibble)
library(DT)

# UI Definition
ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),

  # Application title with styling
  div(
    style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
             padding: 30px; margin: -15px -15px 30px -15px;",
    h1("Cancer Cell Line Gene Analysis",
       style = "color: white; font-weight: 600; margin: 0;"),
    p("Interactive analysis and visualization tool",
      style = "color: rgba(255,255,255,0.9); margin-top: 10px; margin-bottom: 0;")
  ),

  # Main content
  div(
    style = "max-width: 1400px; margin: 0 auto;",

    # Sidebar and main panel
    sidebarLayout(
      # Sidebar Panel
      sidebarPanel(
        width = 3,
        style = "background-color: #f8f9fa; border-radius: 8px; padding: 20px;",

        h4("Analysis Parameters", style = "margin-top: 0; color: #667eea;"),

        # File upload
        div(
          style = "margin-bottom: 25px;",
          fileInput("efx_file",
                    label = div(icon("upload"), " Upload Effect Score File"),
                    accept = c(".csv"),
                    buttonLabel = "Browse",
                    placeholder = "No file selected")
        ),

        # File upload
        div(
          style = "margin-bottom: 25px;",
          fileInput("file",
                    label = div(icon("upload"), " Upload Metadata File"),
                    accept = c(".csv"),
                    buttonLabel = "Browse",
                    placeholder = "No file selected")
        ),

        hr(style = "border-color: #dee2e6;"),

        # Cell line selector - now with multiple selection
        div(
          style = "margin-bottom: 10px;",
          selectizeInput("cellline",
                         label = div(icon("dna"), " Select Cancer Lineage(s)"),
                         choices = NULL,
                         selected = NULL,
                         multiple = TRUE,
                         options = list(
                           placeholder = "Select cancer lineage(s)...",
                           plugins = list("remove_button")
                         ))
        ),

        # Display count of selected cell lines
        div(
          style = "margin-bottom: 20px;",
          uiOutput("cellline_count")
        ),

        # P-value threshold
        div(
          style = "margin-bottom: 20px;",
          numericInput("pvalue",
                       label = div(icon("chart-line"), " P-value Threshold"),
                       value = 0.05,
                       min = 0,
                       max = 1,
                       step = 0.01)
        ),

        # Mean effect score cutoff
        div(
          style = "margin-bottom: 20px;",
          numericInput("mean_effect",
                       label = div(icon("calculator"), " Mean Effect Score Cutoff"),
                       value = 0,
                       max = +Inf,
                       min = -Inf,
                       step = 0.1)
        ),

        # Effect size cutoff
        div(
          style = "margin-bottom: 20px;",
          numericInput("effect_size",
                       label = div(icon("arrows-alt-h"), " Effect Size Cutoff"),
                       value = 0,
                       max = +Inf,
                       min = -Inf,
                       step = 0.1)
        ),

        hr(style = "border-color: #dee2e6;"),

        # Action button
        actionButton("analyze",
                     label = div(icon("play-circle"), " Run Analysis"),
                     class = "btn-primary btn-lg",
                     style = "width: 100%; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                              border: none; font-weight: 500;")
      ),

      # Main Panel
      mainPanel(
        width = 9,

        # Results tabs
        tabsetPanel(
          id = "results_tabs",
          type = "pills",

          # Data Table Tab
          tabPanel(
            title = div(icon("table"), " Results Table"),
            value = "table_tab",
            br(),

            div(
              style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",

              div(
                style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
                h4("Filtered Gene Results", style = "margin: 0; color: #333;"),
                downloadButton("download_table",
                               label = "Download CSV",
                               class = "btn-success",
                               style = "border: none;")
              ),

              hr(style = "margin-top: 10px; margin-bottom: 20px;"),

              DTOutput("results_table")
            )
          ),

          # Volcano Plot Tab
          tabPanel(
            title = div(icon("chart-area"), " Volcano Plot"),
            value = "volcano_tab",
            br(),

            div(
              style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",

              div(
                style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
                h4("Volcano Plot Visualization", style = "margin: 0; color: #333;"),
                downloadButton("download_volcano",
                               label = "Download Plot",
                               class = "btn-success",
                               style = "border: none;")
              ),

              hr(style = "margin-top: 10px; margin-bottom: 20px;"),

              plotOutput("volcano_plot", height = "600px")
            )
          ),

          # Summary Tab
          tabPanel(
            title = div(icon("info-circle"), " Summary"),
            value = "summary_tab",
            br(),

            div(
              style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",

              h4("Analysis Summary", style = "margin-top: 0; color: #333;"),
              hr(),

              uiOutput("summary_stats")
            )
          )
        )
      )
    )
  ),

  # Footer
  div(
    style = "text-align: center; padding: 20px; margin-top: 40px; color: #6c757d; border-top: 1px solid #dee2e6;",
    p("Cancer Cell Line Analysis Tool | Built with Shiny", style = "margin: 0;")
  )
)

# Server Logic
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 500*1024^2)

  # Reactive value to store uploaded data
  efx.data <- reactiveVal(NULL)
  meta.data <- reactiveVal(NULL)

  # Reactive value to store filtered results
  filtered_data <- reactiveVal(NULL)

  # Load and process uploaded file, effect scores
  observeEvent(input$efx_file, {
    req(input$efx_file)

    tryCatch({
      efx <- read.delim(input$efx_file$datapath, sep = ",", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
      # rewrite column names of efx data to be Gene ID only (remove numbers trailing)colnames(efx) <- gsub(" \\([0-9]*\\)", "", colnames(efx))
      colnames(efx) <- gsub(" \\([0-9]*\\)", "", colnames(efx))
      efx.data(efx)

      showNotification("Data loaded successfully!", type = "message", duration = 3)

    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error", duration = 5)
    })
  })

  # Load and process uploaded file, metadata
  observeEvent(input$file, {
    req(input$file)

    tryCatch({
      df <- read.delim(input$file$datapath, sep = ",", row.names = 1, stringsAsFactors = FALSE)
      df <- df[match(rownames(efx.data()), rownames(df)),]
      meta.data(df)

      # Update cell line choices if OncotreeSubtype column exists
      if ("OncotreeSubtype" %in% colnames(df)) {
        cell_lines <- unique(df$OncotreeSubtype)
        updateSelectizeInput(session, "cellline",
                             choices = cell_lines,
                             selected = NULL)
      }

      showNotification("Data loaded successfully!", type = "message", duration = 3)

    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error", duration = 5)
    })
  })

  # Render cell line count display
  output$cellline_count <- renderUI({
    selected <- input$cellline

    if (is.null(selected) || length(selected) == 0) {
      div(
        style = "background: #f0f4f8; padding: 10px; border-radius: 6px; text-align: center;",
        p(style = "margin: 0; color: #6c757d; font-size: 13px;",
          icon("info-circle"), " No lineages selected")
      )
    } else {
      # Count total cell lines across selected lineages
      total_cell_lines <- 0
      if (!is.null(meta.data()) && "OncotreeSubtype" %in% colnames(meta.data())) {
        total_cell_lines <- sum(meta.data()$OncotreeSubtype %in% selected)
      }

      div(
        style = "background: linear-gradient(135deg, #e0e7ff 0%, #ede9fe 100%);
                 padding: 12px; border-radius: 6px; text-align: center;
                 border: 1px solid #c7d2fe;",
        p(style = "margin: 0 0 5px 0; color: #4c1d95; font-weight: 600; font-size: 14px;",
          icon("check-circle"),
          sprintf(" %d lineage%s selected", length(selected), ifelse(length(selected) == 1, "", "s"))),
        p(style = "margin: 0; color: #5b21b6; font-size: 13px;",
          icon("flask"),
          sprintf(" %d total cell line%s", total_cell_lines, ifelse(total_cell_lines == 1, "", "s")))
      )
    }
  })

  # Run analysis when button is clicked
  observeEvent(input$analyze, {
    req(meta.data(), efx.data(), input$cellline)

    if (length(input$cellline) == 0) {
      showNotification("Please select at least one cancer lineage", type = "warning", duration = 3)
      return()
    }
    if (length(rownames(meta.data())) == 0) {
      showNotification("Error: Meta Data not loaded appropriately", type = "error", duration = 3)
      return()
    }
    if (length(rownames(efx.data())) == 0) {
      showNotification("Error: Effects Data not loaded appropriately", type = "error", duration = 3)
      return()
    }

    withProgress(message = 'Running analysis...', value = 0, {

      run_stats <- function(mat, vec, covars = NULL, weights = NULL, target_type = 'Gene') {
          require(limma)
          require(plyr)
          require(magrittr)
          require(dplyr)
          require(tibble)
          
          udata <- which(!is.na(vec))
          if (!is.numeric(vec)) {
            pred <- factor(vec[udata])
            stopifnot(length(levels(pred)) == 2) #only two group comparisons implemented so far
            n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],,drop=F]))
            n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],,drop=F]))
            min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
          } else {
            pred <- vec[udata]
            min_samples <- colSums(!is.na(mat[udata,]))
          }
          #there must be more than one unique value of the independent variable
          if (length(unique(pred)) <= 1) {
            return(NULL)
          }
          #if using covariates add them as additional predictors to the model
          if (!is.null(covars)) {
            if (!is.data.frame(covars)) {
              covars <- data.frame(covars)
            }
            combined <- covars[udata,, drop = FALSE]
            combined[['pred']] <- pred
            form <- as.formula(paste('~', paste0(colnames(combined), collapse = ' + ')))
            design <- model.matrix(form, combined)
            design <- design[, colSums(design) != 0, drop = FALSE]
          } else {
            design <- model.matrix(~pred)
          }
          if (!is.null(weights)) {
            if (is.matrix(weights)) {
              weights <- t(weights[udata,])
            } else{
              weights <- weights[udata]
            }
          }
          fit <- limma::lmFit(t(mat[udata,]), design, weights = weights)
          fit <- limma::eBayes(fit)
          targ_coef <- grep('pred', colnames(design), value = TRUE)
          res <- limma::topTable(fit, coef = targ_coef, number = Inf)
          
          if (colnames(res)[1] == 'ID') {
            colnames(res)[1] <- target_type
          } else {
            res %<>% rownames_to_column(var = target_type)
          }
          res$min_samples <- min_samples[res[[target_type]]]
          
          two_to_one_sided <- function(two_sided_p, stat, test_dir) {
            #helper function for converting two-sided p-values to one-sided p-values
            one_sided_p <- two_sided_p / 2
            if (test_dir == 'right') {
              one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 0]
            } else {
              one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 0]
            }
            return(one_sided_p)
          }
          res %<>% 
            set_colnames(
              revalue(
                colnames(.), 
                c('logFC' = 'EffectSize', 'AveExpr' = 'Avg', 't' = 't_stat', 'B' = 'log_odds', 
                'P.Value' = 'p.value', 'adj.P.Val' = 'q.value', 'min_samples' = 'min_samples')
                )
              ) %>% 
                na.omit()
          res %<>% 
            dplyr::mutate(p.left = two_to_one_sided(p.value, EffectSize, 'left'),
                                    p.right = two_to_one_sided(p.value, EffectSize, 'right'),
                                    q.left = p.adjust(p.left, method = 'BH'),
                                    q.right = p.adjust(p.right, method = 'BH'))
          return(res)
        }

      # Create binary indicator: 1 if cell line matches selected lineages, 0 otherwise
      ind.var <- as.numeric(meta.data()$OncotreeSubtype %in% input$cellline)

      # Run limma analysis comparing selected lineages vs others
      res1 <- run_stats(
        mat = efx.data() %>% as.matrix(),
        vec = ind.var
      )

      # Check if limma returned results
      if (length(rownames(res1)) == 0) {
        showNotification("Analysis failed: insufficient variation in data", type = "error", duration = 10)
        return()
      }

      # Get the effect scores for only the selected lineages to calculate mean
      effects_subset <- efx.data() %>%
        subset(rownames(efx.data()) %in% rownames(meta.data())[which(meta.data()$OncotreeSubtype %in% input$cellline)])

      # Calculate mean effect scores across selected cell lines
      mean_effects <- effects_subset %>%
        map_dbl(~mean(as.numeric(.))) %>%
        enframe(name = "Gene", value = "Mean_Effect")

      # Join and filter based on user inputs
      results <- res1 %>%
        left_join(mean_effects, by = "Gene") %>%
        filter(q.value < input$pvalue,
               Mean_Effect < input$mean_effect,
               EffectSize < input$effect_size)

      incProgress(1, detail = "Analysis complete")
      filtered_data(results)

      showNotification(
        paste("Analysis complete!", nrow(results), "genes meet the criteria"),
        type = "message",
        duration = 4
      )
    })
  })

  # Render results table
  output$results_table <- renderDT({
    req(filtered_data())

    datatable(
      filtered_data(),
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        searchHighlight = TRUE,
        ordering = TRUE
      ),
      class = 'cell-border stripe hover',
      rownames = FALSE,
      filter = 'top'
    )
  })

  # Generate volcano plot
  output$volcano_plot <- renderPlot({
    req(filtered_data())

    # Create volcano plot using actual results data
    plot_data <- filtered_data() %>%
      mutate(neg_log10_pval = -log10(q.value))

    ggplot(plot_data, aes(x = EffectSize, y = neg_log10_pval)) +
      geom_point(aes(color = Mean_Effect), alpha = 0.7, size = 3) +
      geom_hline(yintercept = -log10(input$pvalue), linetype = "dashed", color = "red", alpha = 0.7) +
      geom_vline(xintercept = input$effect_size, linetype = "dashed", color = "blue", alpha = 0.7) +
      scale_color_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                            midpoint = median(plot_data$Mean_Effect, na.rm = TRUE),
                            name = "Mean Effect") +
      labs(
        title = "Volcano Plot",
        subtitle = paste("Cell Lineages:", paste(input$cellline, collapse = ", ")),
        x = "Effect Size",
        y = "-Log10(P-value)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 18, color = "#333"),
        plot.subtitle = element_text(color = "#666", margin = margin(b = 15)),
        panel.grid.major = element_line(color = "#e0e0e0"),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold", color = "#555"),
        legend.position = "right"
      ) +
      coord_flip()
  })

  # Generate summary statistics
  output$summary_stats <- renderUI({
    req(filtered_data())

    # ============================================
    # PLACEHOLDER: Add your summary statistics here
    # ============================================
    # Calculate relevant statistics from your filtered data

    div(
      div(
        class = "row",
        div(
          class = "col-md-6",
          div(
            style = "background: #f0f7ff; padding: 15px; border-radius: 6px; margin-bottom: 15px;",
            h5(icon("dna"), " Selected Cell Lineages", style = "color: #667eea; margin-top: 0;"),
            p(strong(paste(input$cellline, collapse = ", ")), style = "font-size: 14px; margin: 0;"),
            p(em(paste(length(input$cellline), "lineage(s) selected")),
              style = "font-size: 12px; color: #6c757d; margin: 5px 0 0 0;")
          )
        ),
        div(
          class = "col-md-6",
          div(
            style = "background: #f0fff4; padding: 15px; border-radius: 6px; margin-bottom: 15px;",
            h5(icon("list"), " Total Genes Passing Filters", style = "color: #48bb78; margin-top: 0;"),
            p(strong(nrow(filtered_data())), style = "font-size: 16px; margin: 0;")
          )
        )
      ),

      hr(),

      h5("Filter Criteria Applied:", style = "margin-top: 20px;"),
      tags$ul(
        style = "line-height: 2;",
        tags$li(strong("P-value threshold: "), input$pvalue),
        tags$li(strong("Mean effect score cutoff: "), input$mean_effect),
        tags$li(strong("Effect size cutoff: "), input$effect_size)
      ),

      div(
        style = "background: #fffbeb; border-left: 4px solid #f59e0b; padding: 15px; margin-top: 20px; border-radius: 4px;",
        p(icon("info-circle"), strong(" Note: "),
          "Additional summary statistics will appear here based on your analysis results.",
          style = "margin: 0; color: #92400e;")
      )
    )
  })

  # Download filtered table
  output$download_table <- downloadHandler(
    filename = function() {
      paste0("filtered_genes_", paste(input$cellline, collapse = "_"), "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(filtered_data())
      write.csv(filtered_data(), file, row.names = FALSE)
    }
  )

  # Download volcano plot
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0("volcano_plot_", paste(input$cellline, collapse = "_"), "_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(filtered_data())

      # ============================================
      # PLACEHOLDER: Save your actual volcano plot here
      # ============================================

      ggsave(
        file,
        plot = last_plot(),
        width = 12,
        height = 8,
        dpi = 300,
        bg = "white"
      )
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
