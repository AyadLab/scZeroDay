library(shiny)
library(tidyverse)
library(plyr)
library(dplyr)
library(tibble)
library(DT)
library(plotly)
library(enrichR)
library(openxlsx)

# UI Definition
ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),

  # Application title with styling
  div(
    style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
             padding: 30px; margin: -15px -15px 30px -15px;",
    h1("scZeroDay",
       style = "color: white; font-weight: 600; margin: 0;"),
    p("Interactive gene essentiality analysis and visualization tool",
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
                    label = div(icon("upload"), " (1) Upload Effect Score File"),
                    accept = c(".csv"),
                    buttonLabel = "Browse",
                    placeholder = "No file selected")
        ),

        # File upload
        div(
          style = "margin-bottom: 25px;",
          fileInput("file",
                    label = div(icon("upload"), " (2) Upload Metadata File"),
                    accept = c(".csv"),
                    buttonLabel = "Browse",
                    placeholder = "No file selected")
        ),

        # Lineage column selector
        div(
          style = "margin-bottom: 25px;",
          selectInput("lineage_column",
                      label = div(icon("columns"), " (3) Select Lineage Column"),
                      choices = NULL,
                      selected = NULL)
        ),

        hr(style = "border-color: #dee2e6;"),

        # Cell line selector - now with multiple selection
        div(
          style = "margin-bottom: 10px;",
          selectizeInput("cellline",
                         label = div(icon("dna"), " (4) Select Lineage(s) of Interest"),
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
                       label = div(icon("chart-line"), " (5) P-value Threshold"),
                       value = 0.05,
                       min = 0,
                       max = 1,
                       step = 0.01)
        ),

        # Effect size cutoff
        div(
          style = "margin-bottom: 20px;",
          numericInput("effect_size",
                       label = div(icon("arrows-alt-h"), " (6) Effect Size Cutoff"),
                       value = 0,
                       max = +Inf,
                       min = -Inf,
                       step = 0.1)
        ),

        # Mean effect score cutoff
        div(
          style = "margin-bottom: 20px;",
          numericInput("mean_effect",
                       label = div(icon("calculator"), " (7) Mean Effect Score Cutoff"),
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

              plotlyOutput("volcano_plot", height = "600px")
            )
          ),

          # Pathway Enrichment Tab
          tabPanel(
            title = div(icon("project-diagram"), " Pathway Enrichment"),
            value = "enrichment_tab",
            br(),

            div(
              style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",

              div(
                style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
                h4("Pathway Enrichment Analysis", style = "margin: 0; color: #333;"),
                downloadButton("download_enrichment",
                               label = "Download Results",
                               class = "btn-success",
                               style = "border: none;")
              ),

              hr(style = "margin-top: 10px; margin-bottom: 20px;"),

              p("Enrichment analysis of filtered genes using multiple pathway databases:",
                style = "color: #666; margin-bottom: 20px;"),

              tabsetPanel(
                id = "enrichment_tabs",
                type = "tabs",

                tabPanel(
                  title = "KEGG",
                  br(),
                  DTOutput("enrichr_kegg")
                ),

                tabPanel(
                  title = "WikiPathways",
                  br(),
                  DTOutput("enrichr_wiki")
                ),

                tabPanel(
                  title = "Reactome",
                  br(),
                  DTOutput("enrichr_reactome")
                ),

                tabPanel(
                  title = "MSigDB Hallmark",
                  br(),
                  DTOutput("enrichr_msigdb")
                )
              )
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
    p("scZeroDay | Built with Shiny", style = "margin: 0;")
  )
)

# Server Logic
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 500*1024^2)

  # Reactive value to store uploaded data
  efx.data <- reactiveVal(NULL)
  meta.data <- reactiveVal(NULL)

  # Reactive value to store filtered results
  raw <- reactiveVal(NULL)
  filtered_data <- reactiveVal(NULL)

  # Reactive value to store enrichment results
  enrichment_results <- reactiveVal(NULL)

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

      # Update lineage column choices with all column names
      updateSelectInput(session, "lineage_column",
                        choices = colnames(df),
                        selected = if ("OncotreeSubtype" %in% colnames(df)) "OncotreeSubtype" else colnames(df)[1])

      showNotification("Metadata loaded successfully!", type = "message", duration = 3)

    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error", duration = 5)
    })
  })

  # Update cell line choices when lineage column is selected
  observeEvent(input$lineage_column, {
    req(meta.data(), input$lineage_column)

    tryCatch({
      df <- meta.data()
      if (input$lineage_column %in% colnames(df)) {
        cell_lines <- unique(df[[input$lineage_column]])
        updateSelectizeInput(session, "cellline",
                             choices = cell_lines,
                             selected = NULL)
      }
    }, error = function(e) {
      showNotification(paste("Error updating cell line choices:", e$message), type = "error", duration = 5)
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
      if (!is.null(meta.data()) && !is.null(input$lineage_column) && input$lineage_column %in% colnames(meta.data())) {
        total_cell_lines <- sum(meta.data()[[input$lineage_column]] %in% selected)
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
    req(meta.data(), efx.data(), input$cellline, input$lineage_column)

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
    if (!input$lineage_column %in% colnames(meta.data())) {
      showNotification("Error: Selected lineage column not found in metadata", type = "error", duration = 3)
      return()
    }

    # Check cell line count and warn if fewer than 10
    total_cell_lines <- sum(meta.data()[[input$lineage_column]] %in% input$cellline)
    if (total_cell_lines < 10) {
      showNotification(
        paste0("Warning: Only ", total_cell_lines, " cell line(s) selected. ",
               "Fewer than 10 cell lines may impact statistical validity."),
        type = "warning",
        duration = 8
      )
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
      ind.var <- as.numeric(meta.data()[[input$lineage_column]] %in% input$cellline)

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
        subset(rownames(efx.data()) %in% rownames(meta.data())[which(meta.data()[[input$lineage_column]] %in% input$cellline)])

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

      res2 <- res1 %>%
        left_join(mean_effects, by = "Gene")

      incProgress(0.7, detail = "Running pathway enrichment...")

      # Run enrichment analysis if there are filtered genes
      enrichment_data <- NULL
      if (nrow(results) > 0) {
        tryCatch({
          gene_list <- results$Gene

          # Run enrichR on multiple databases
          dbs <- c("KEGG_2026", "WikiPathway_2024_Human", "Reactome_Pathways_2024", "MSigDB_Hallmark_2020") ### Update database names to match enrichR's current offerings
          enrichment_data <- enrichr(gene_list, dbs)

          # Rename the list elements for clearer display
          names(enrichment_data) <- c("KEGG_2026", "WikiPathway_2024", "Reactome_2024", "MSigDB_Hallmark_2020")

        }, error = function(e) {
          showNotification(
            paste("Warning: Enrichment analysis failed:", e$message),
            type = "warning",
            duration = 5
          )
        })
      }

      incProgress(1, detail = "Analysis complete")
      raw(res2)
      filtered_data(results)
      enrichment_results(enrichment_data)

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

  # Generate interactive volcano plot
  output$volcano_plot <- renderPlotly({
    req(filtered_data(), raw())

    # Create volcano plot using actual results data
    plot_data <- raw() %>%
      mutate(neg_log10_pval = -log10(q.value),
             # Create custom hover text
             hover_text = paste0(
               "<b>Gene:</b> ", Gene, "<br>",
               "<b>Effect Size:</b> ", round(EffectSize, 4), "<br>",
               "<b>Mean Effect:</b> ", round(Mean_Effect, 4), "<br>",
               "<b>q-value:</b> ", formatC(q.value, format = "e", digits = 3)
             ))

    # Identify filtered/highlighted points
    filtered_genes <- filtered_data()$Gene

    p <- ggplot(plot_data, aes(x = EffectSize, y = neg_log10_pval, text = hover_text)) +
      geom_point(aes(color = Mean_Effect), alpha = 0.7, size = 3) +
      geom_hline(yintercept = -log10(input$pvalue), linetype = "dashed", color = "red", alpha = 0.7) +
      geom_vline(xintercept = input$effect_size, linetype = "dashed", color = "blue", alpha = 0.7) +
      scale_color_gradient2(low = "#4400ff", mid = "#fafafa", high = "#ff2c02",
                            midpoint = median(plot_data$Mean_Effect, na.rm = TRUE),
                            name = "Mean Effect") +
      labs(
        title = "Volcano Plot",
        subtitle = paste("Cell Lineages:", paste(input$cellline, collapse = ", ")),
        x = "Effect Size",
        y = "-Log10(q-value)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 18, color = "#333"),
        plot.subtitle = element_text(color = "#666", margin = margin(b = 15)),
        panel.grid.major = element_line(color = "#e0e0e0"),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold", color = "#555"),
        legend.position = "right"
      )

    # Convert to plotly with custom tooltip
    ggplotly(p, tooltip = "text") %>%
      layout(
        hoverlabel = list(
          bgcolor = "white",
          font = list(size = 12, color = "black"),
          bordercolor = "#667eea"
        )
      )
  })

  # Generate summary statistics
  output$summary_stats <- renderUI({
    req(filtered_data(), raw())

    # Calculate summary statistics
    filtered <- filtered_data()
    all_data <- raw()

    # Total genes analyzed
    total_genes <- nrow(all_data)
    filtered_genes <- nrow(filtered)

    # Effect size statistics for filtered genes
    if (filtered_genes > 0) {
      mean_effect_size <- mean(filtered$EffectSize, na.rm = TRUE)
      median_effect_size <- median(filtered$EffectSize, na.rm = TRUE)
      min_effect_size <- min(filtered$EffectSize, na.rm = TRUE)
      max_effect_size <- max(filtered$EffectSize, na.rm = TRUE)

      mean_mean_effect <- mean(filtered$Mean_Effect, na.rm = TRUE)
      median_mean_effect <- median(filtered$Mean_Effect, na.rm = TRUE)

      min_qvalue <- min(filtered$q.value, na.rm = TRUE)
      max_qvalue <- max(filtered$q.value, na.rm = TRUE)

      # Top 5 genes by effect size (most negative)
      top_genes <- filtered %>%
        arrange(EffectSize) %>%
        head(5)

      # Top 5 genes by mean effect score (most negative)
      top_genes_mean_effect <- filtered %>%
        arrange(Mean_Effect) %>%
        head(5)
    } else {
      mean_effect_size <- median_effect_size <- min_effect_size <- max_effect_size <- NA
      mean_mean_effect <- median_mean_effect <- NA
      min_qvalue <- max_qvalue <- NA
      top_genes <- NULL
      top_genes_mean_effect <- NULL
    }

    # Count cell lines in selected lineages
    total_cell_lines <- 0
    if (!is.null(meta.data()) && !is.null(input$lineage_column) && input$lineage_column %in% colnames(meta.data())) {
      total_cell_lines <- sum(meta.data()[[input$lineage_column]] %in% input$cellline)
    }

    div(
      # Row 1: Overview cards
      div(
        class = "row",
        div(
          class = "col-md-4",
          div(
            style = "background: #f0f7ff; padding: 15px; border-radius: 6px; margin-bottom: 15px;",
            h5(icon("dna"), " Selected Cell Lineages", style = "color: #667eea; margin-top: 0;"),
            p(strong(paste(input$cellline, collapse = ", ")), style = "font-size: 14px; margin: 0;"),
            p(em(paste(length(input$cellline), "lineage(s),", total_cell_lines, "cell line(s)")),
              style = "font-size: 12px; color: #6c757d; margin: 5px 0 0 0;")
          )
        ),
        div(
          class = "col-md-4",
          div(
            style = "background: #f0fff4; padding: 15px; border-radius: 6px; margin-bottom: 15px;",
            h5(icon("list"), " Genes Passing Filters", style = "color: #48bb78; margin-top: 0;"),
            p(strong(filtered_genes), " of ", total_genes, " total genes",
              style = "font-size: 16px; margin: 0;"),
            p(em(paste0(round(100 * filtered_genes / total_genes, 2), "% of total")),
              style = "font-size: 12px; color: #6c757d; margin: 5px 0 0 0;")
          )
        ),
        div(
          class = "col-md-4",
          div(
            style = "background: #fef3c7; padding: 15px; border-radius: 6px; margin-bottom: 15px;",
            h5(icon("chart-bar"), " q-value Range", style = "color: #d97706; margin-top: 0;"),
            if (filtered_genes > 0) {
              tagList(
                p(strong("Min: "), formatC(min_qvalue, format = "e", digits = 2),
                  style = "font-size: 14px; margin: 0;"),
                p(strong("Max: "), formatC(max_qvalue, format = "e", digits = 2),
                  style = "font-size: 14px; margin: 5px 0 0 0;")
              )
            } else {
              p("No data", style = "font-size: 14px; margin: 0; color: #6c757d;")
            }
          )
        )
      ),

      hr(),

      # Row 2: Effect statistics
      div(
        class = "row",
        div(
          class = "col-md-6",
          div(
            style = "background: #fdf2f8; padding: 15px; border-radius: 6px; margin-bottom: 15px;",
            h5(icon("arrows-alt-h"), " Effect Size Statistics (Filtered)", style = "color: #be185d; margin-top: 0;"),
            if (filtered_genes > 0) {
              tags$table(
                style = "width: 100%; font-size: 14px;",
                tags$tr(tags$td(strong("Mean:")), tags$td(round(mean_effect_size, 4))),
                tags$tr(tags$td(strong("Median:")), tags$td(round(median_effect_size, 4))),
                tags$tr(tags$td(strong("Min:")), tags$td(round(min_effect_size, 4))),
                tags$tr(tags$td(strong("Max:")), tags$td(round(max_effect_size, 4)))
              )
            } else {
              p("No filtered genes to display", style = "margin: 0; color: #6c757d;")
            }
          )
        ),
        div(
          class = "col-md-6",
          div(
            style = "background: #ecfdf5; padding: 15px; border-radius: 6px; margin-bottom: 15px;",
            h5(icon("calculator"), " Mean Effect Statistics (Filtered)", style = "color: #059669; margin-top: 0;"),
            if (filtered_genes > 0) {
              tags$table(
                style = "width: 100%; font-size: 14px;",
                tags$tr(tags$td(strong("Mean:")), tags$td(round(mean_mean_effect, 4))),
                tags$tr(tags$td(strong("Median:")), tags$td(round(median_mean_effect, 4)))
              )
            } else {
              p("No filtered genes to display", style = "margin: 0; color: #6c757d;")
            }
          )
        )
      ),

      hr(),

      # Row 3: Top genes table
      if (filtered_genes > 0 && !is.null(top_genes) && nrow(top_genes) > 0) {
        div(
          style = "background: #f8fafc; padding: 15px; border-radius: 6px; margin-bottom: 15px;",
          h5(icon("trophy"), " Top 5 Genes by Effect Size", style = "color: #475569; margin-top: 0;"),
          tags$table(
            class = "table table-striped table-sm",
            style = "font-size: 13px; margin-bottom: 0;",
            tags$thead(
              tags$tr(
                tags$th("Gene"),
                tags$th("Effect Size"),
                tags$th("Mean Effect"),
                tags$th("q-value")
              )
            ),
            tags$tbody(
              lapply(1:nrow(top_genes), function(i) {
                tags$tr(
                  tags$td(strong(top_genes$Gene[i])),
                  tags$td(round(top_genes$EffectSize[i], 4)),
                  tags$td(round(top_genes$Mean_Effect[i], 4)),
                  tags$td(formatC(top_genes$q.value[i], format = "e", digits = 2))
                )
              })
            )
          )
        )
      },

      # Row 4: Top genes table, effect score
      if (filtered_genes > 0 && !is.null(top_genes_mean_effect) && nrow(top_genes_mean_effect) > 0) {
        div(
          style = "background: #f8fafc; padding: 15px; border-radius: 6px; margin-bottom: 15px;",
          h5(icon("trophy"), " Top 5 Genes by Mean Effect Score", style = "color: #475569; margin-top: 0;"),
          tags$table(
            class = "table table-striped table-sm",
            style = "font-size: 13px; margin-bottom: 0;",
            tags$thead(
              tags$tr(
                tags$th("Gene"),
                tags$th("Effect Size"),
                tags$th("Mean Effect"),
                tags$th("q-value")
              )
            ),
            tags$tbody(
              lapply(1:nrow(top_genes_mean_effect), function(i) {
                tags$tr(
                  tags$td(strong(top_genes_mean_effect$Gene[i])),
                  tags$td(round(top_genes_mean_effect$EffectSize[i], 4)),
                  tags$td(round(top_genes_mean_effect$Mean_Effect[i], 4)),
                  tags$td(formatC(top_genes_mean_effect$q.value[i], format = "e", digits = 2))
                )
              })
            )
          )
        )
      },

      hr(),

      h5("Filter Criteria Applied:", style = "margin-top: 20px;"),
      tags$ul(
        style = "line-height: 2;",
        tags$li(strong("q-value threshold: "), input$pvalue),
        tags$li(strong("Mean effect score cutoff: "), input$mean_effect),
        tags$li(strong("Effect size cutoff: "), input$effect_size)
      )
    )
  })

  # Render enrichment tables
  output$enrichr_kegg <- renderDT({
    req(enrichment_results())
    enrich_data <- enrichment_results()

    if (!is.null(enrich_data) && "KEGG_2026" %in% names(enrich_data)) {
      df <- enrich_data$KEGG_2026 %>%
        filter(Adjusted.P.value < 0.05) %>%
        select(Term, Overlap, P.value, Adjusted.P.value, Genes) %>%
        arrange(P.value)

      datatable(
        df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          searchHighlight = TRUE,
          ordering = TRUE
        ),
        class = 'cell-border stripe hover',
        rownames = FALSE,
        filter = 'top'
      ) %>%
        formatSignif(columns = c('P.value', 'Adjusted.P.value'), digits = 3)
    } else {
      datatable(data.frame(Message = "No enrichment results available"))
    }
  })

  output$enrichr_wiki <- renderDT({
    req(enrichment_results())
    enrich_data <- enrichment_results()

    if (!is.null(enrich_data) && "WikiPathways_2024" %in% names(enrich_data)) {
      df <- enrich_data$WikiPathways_2024 %>%
        filter(Adjusted.P.value < 0.05) %>%
        select(Term, Overlap, P.value, Adjusted.P.value, Genes) %>%
        arrange(P.value)

      datatable(
        df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          searchHighlight = TRUE,
          ordering = TRUE
        ),
        class = 'cell-border stripe hover',
        rownames = FALSE,
        filter = 'top'
      ) %>%
        formatSignif(columns = c('P.value', 'Adjusted.P.value'), digits = 3)
    } else {
      datatable(data.frame(Message = "No enrichment results available"))
    }
  })

  output$enrichr_reactome <- renderDT({
    req(enrichment_results())
    enrich_data <- enrichment_results()

    if (!is.null(enrich_data) && "Reactome_2024" %in% names(enrich_data)) {
      df <- enrich_data$Reactome_2024 %>%
        filter(Adjusted.P.value < 0.05) %>%
        select(Term, Overlap, P.value, Adjusted.P.value, Genes) %>%
        arrange(P.value)

      datatable(
        df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          searchHighlight = TRUE,
          ordering = TRUE
        ),
        class = 'cell-border stripe hover',
        rownames = FALSE,
        filter = 'top'
      ) %>%
        formatSignif(columns = c('P.value', 'Adjusted.P.value'), digits = 3)
    } else {
      datatable(data.frame(Message = "No enrichment results available"))
    }
  })

  output$enrichr_msigdb <- renderDT({
    req(enrichment_results())
    enrich_data <- enrichment_results()

    if (!is.null(enrich_data) && "MSigDB_Hallmark_2020" %in% names(enrich_data)) {
      df <- enrich_data$MSigDB_Hallmark_2020 %>%
        filter(Adjusted.P.value < 0.05) %>%
        select(Term, Overlap, P.value, Adjusted.P.value, Genes) %>%
        arrange(P.value)

      datatable(
        df,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          searchHighlight = TRUE,
          ordering = TRUE
        ),
        class = 'cell-border stripe hover',
        rownames = FALSE,
        filter = 'top'
      ) %>%
        formatSignif(columns = c('P.value', 'Adjusted.P.value'), digits = 3)
    } else {
      datatable(data.frame(Message = "No enrichment results available"))
    }
  })

  # Download enrichment results
  output$download_enrichment <- downloadHandler(
    filename = function() {
      paste0("enrichment_results_", paste(input$cellline, collapse = "_"), "_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(enrichment_results())

      # Create a workbook and add each enrichment database as a sheet
      library(openxlsx)
      wb <- createWorkbook()

      enrich_data <- enrichment_results()
      if (!is.null(enrich_data)) {
        for (db_name in names(enrich_data)) {
          addWorksheet(wb, db_name)
          writeData(wb, db_name, enrich_data[[db_name]])
        }
      }

      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

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
      req(filtered_data(), raw())

      # Recreate the volcano plot for saving
      plot_data <- raw() %>%
        mutate(neg_log10_pval = -log10(q.value))

      # Identify filtered/highlighted points
      filtered_genes <- filtered_data()$Gene

      p <- ggplot(plot_data, aes(x = EffectSize, y = neg_log10_pval)) +
        geom_point(aes(color = Mean_Effect), alpha = 0.7, size = 3) +
        geom_hline(yintercept = -log10(input$pvalue), linetype = "dashed", color = "red", alpha = 0.7) +
        geom_vline(xintercept = input$effect_size, linetype = "dashed", color = "blue", alpha = 0.7) +
        geom_label(data = plot_data %>% filter(Gene %in% filtered_genes),
                   aes(label = Gene),
                   size = 3, color = "black", fill = "white",
                   alpha = 0.8, linewidth = 0.2, position = "nudge") +
        scale_color_gradient2(low = "#1974d0", mid = "#f7f7f78e", high = "#e30a23",
                              midpoint = median(plot_data$Mean_Effect, na.rm = TRUE),
                              name = "Mean Effect") +
        labs(
          title = "Volcano Plot",
          subtitle = paste("Cell Lineages:", paste(input$cellline, collapse = ", ")),
          x = "Effect Size",
          y = "-Log10(q-value)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", size = 18, color = "#333"),
          plot.subtitle = element_text(color = "#666", margin = margin(b = 15)),
          panel.grid.major = element_line(color = "#e0e0e0"),
          panel.grid.minor = element_blank(),
          axis.title = element_text(face = "bold", color = "#555"),
          legend.position = "right"
        )

      ggsave(
        file,
        plot = p,
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
