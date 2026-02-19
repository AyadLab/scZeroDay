library(shiny)
library(tidyverse)
library(plyr)
library(dplyr)
library(tibble)
library(DT)
library(plotly)
library(enrichR)
library(openxlsx)
library(httr)
library(jsonlite)

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
                  title = "WikiPathway",
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

          # ChEA3 Transcription Factor Analysis Tab
          tabPanel(
            title = div(icon("dna"), " ChEA3 TF Analysis"),
            value = "chea3_tab",
            br(),

            div(
              style = "background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",

              div(
                style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
                h4("ChEA3 Transcription Factor Enrichment Analysis", style = "margin: 0; color: #333;"),
                downloadButton("download_chea3",
                               label = "Download Results",
                               class = "btn-success",
                               style = "border: none;")
              ),

              hr(style = "margin-top: 10px; margin-bottom: 20px;"),

              p("Identify transcription factors that may regulate your gene set using ",
                tags$a(href = "https://maayanlab.cloud/chea3/", target = "_blank", "ChEA3"),
                " (ChIP-X Enrichment Analysis 3).",
                style = "color: #666; margin-bottom: 20px;"),

              # Gene source selection
              div(
                style = "background: #f8f9fa; padding: 15px; border-radius: 6px; margin-bottom: 20px;",

                h5("Select Gene Source", style = "margin-top: 0; color: #667eea;"),

                tags$style("#chea3_gene_source .radio-inline { margin-right: 25px; }"),
                radioButtons("chea3_gene_source",
                             label = NULL,
                             choices = c("Filtered Gene Set" = "filtered",
                                         "EnrichR Pathway Genes" = "enrichr"),
                             selected = "filtered",
                             inline = TRUE),

                # Conditional panel for EnrichR pathway selection
                conditionalPanel(
                  condition = "input.chea3_gene_source == 'enrichr'",

                  div(
                    style = "margin-top: 15px;",

                    selectInput("chea3_enrichr_db",
                                label = "Select Pathway Database:",
                                choices = c("KEGG" = "KEGG_2026",
                                            "WikiPathway" = "WikiPathways_2024",
                                            "Reactome" = "Reactome_2024",
                                            "MSigDB Hallmark" = "MSigDB_Hallmark_2020"),
                                selected = "KEGG_2026"),

                    selectInput("chea3_pathway",
                                label = "Select Pathway:",
                                choices = NULL,
                                selected = NULL)
                  )
                ),

                # Display selected genes info
                div(
                  style = "margin-top: 15px;",
                  uiOutput("chea3_gene_info")
                ),

                # Run ChEA3 button
                div(
                  style = "margin-top: 15px;",
                  actionButton("run_chea3",
                               label = div(icon("search"), " Run ChEA3 Analysis"),
                               class = "btn-primary",
                               style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); border: none;")
                )
              ),

              # Results section
              div(
                style = "margin-top: 20px;",

                tabsetPanel(
                  id = "chea3_results_tabs",
                  type = "tabs",

                  tabPanel(
                    title = "Integrated Results (Mean Rank)",
                    br(),
                    DTOutput("chea3_integrated")
                  ),

                  tabPanel(
                    title = "ENCODE",
                    br(),
                    DTOutput("chea3_encode")
                  ),

                  tabPanel(
                    title = "ReMap",
                    br(),
                    DTOutput("chea3_remap")
                  ),

                  tabPanel(
                    title = "ARCHS4",
                    br(),
                    DTOutput("chea3_archs4")
                  ),

                  tabPanel(
                    title = "GTEx",
                    br(),
                    DTOutput("chea3_gtex")
                  ),

                  tabPanel(
                    title = "Literature",
                    br(),
                    DTOutput("chea3_literature")
                  ),

                  tabPanel(
                    title = "Enrichr",
                    br(),
                    DTOutput("chea3_enrichr")
                  )
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

  # Reactive value to store ChEA3 results
  chea3_results <- reactiveVal(NULL)

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
          dbs <- c("KEGG_2026", "WikiPathways_2024_Human", "Reactome_Pathways_2024", "MSigDB_Hallmark_2020") ### Update database names to match enrichR's current offerings
          enrichment_data <- enrichr(gene_list, dbs)

          # Rename the list elements for clearer display
          names(enrichment_data) <- c("KEGG_2026", "WikiPathways_2024", "Reactome_2024", "MSigDB_Hallmark_2020")

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

  # ChEA3 Tab Logic

  # Update pathway choices when database is selected
  observeEvent(list(input$chea3_enrichr_db, enrichment_results()), {
    req(enrichment_results(), input$chea3_enrichr_db)

    enrich_data <- enrichment_results()
    db_name <- input$chea3_enrichr_db

    if (!is.null(enrich_data) && db_name %in% names(enrich_data)) {
      # Get significant pathways (Adjusted.P.value < 0.05)
      pathways <- enrich_data[[db_name]] %>%
        filter(Adjusted.P.value < 0.05) %>%
        arrange(P.value) %>%
        pull(Term)

      if (length(pathways) > 0) {
        updateSelectInput(session, "chea3_pathway",
                          choices = pathways,
                          selected = pathways[1])
      } else {
        updateSelectInput(session, "chea3_pathway",
                          choices = c("No significant pathways found" = ""),
                          selected = "")
      }
    } else {
      updateSelectInput(session, "chea3_pathway",
                        choices = c("Run analysis first" = ""),
                        selected = "")
    }
  })

  # Reactive to get genes for ChEA3 based on source selection
  chea3_genes <- reactive({
    if (input$chea3_gene_source == "filtered") {
      req(filtered_data())
      return(filtered_data()$Gene)
    } else {
      req(enrichment_results(), input$chea3_enrichr_db, input$chea3_pathway)

      enrich_data <- enrichment_results()
      db_name <- input$chea3_enrichr_db

      if (!is.null(enrich_data) && db_name %in% names(enrich_data) && input$chea3_pathway != "") {
        # Get genes from the selected pathway
        pathway_row <- enrich_data[[db_name]] %>%
          filter(Term == input$chea3_pathway)

        if (nrow(pathway_row) > 0) {
          genes <- pathway_row$Genes[1]
          # Split by semicolon and clean up
          gene_list <- unlist(strsplit(genes, ";"))
          gene_list <- trimws(gene_list)
          return(gene_list[gene_list != ""])
        }
      }
      return(NULL)
    }
  })

  # Render gene info display
  output$chea3_gene_info <- renderUI({
    genes <- chea3_genes()

    if (is.null(genes) || length(genes) == 0) {
      div(
        style = "background: #fef3c7; padding: 10px; border-radius: 6px;",
        p(style = "margin: 0; color: #d97706;",
          icon("exclamation-triangle"), " No genes available. Run the main analysis first or select a valid pathway.")
      )
    } else {
      div(
        style = "background: #e0e7ff; padding: 10px; border-radius: 6px;",
        p(style = "margin: 0; color: #4c1d95; font-weight: 600;",
          icon("check-circle"),
          sprintf(" %d genes ready for ChEA3 analysis", length(genes))),
        if (length(genes) <= 10) {
          p(style = "margin: 5px 0 0 0; color: #5b21b6; font-size: 12px;",
            paste(genes, collapse = ", "))
        } else {
          p(style = "margin: 5px 0 0 0; color: #5b21b6; font-size: 12px;",
            paste(c(head(genes, 10), "..."), collapse = ", "))
        }
      )
    }
  })

  # Run ChEA3 analysis
  observeEvent(input$run_chea3, {
    genes <- chea3_genes()

    if (is.null(genes) || length(genes) == 0) {
      showNotification("No genes available for ChEA3 analysis", type = "error", duration = 3)
      return()
    }

    if (length(genes) < 3) {
      showNotification("ChEA3 requires at least 3 genes", type = "error", duration = 3)
      return()
    }

    withProgress(message = 'Running ChEA3 analysis...', value = 0, {
      tryCatch({
        # Prepare the request body
        body <- list(
          query_name = "shiny_query",
          gene_set = genes
        )

        incProgress(0.3, detail = "Sending request to ChEA3...")

        # Make the API call
        response <- POST(
          url = "https://maayanlab.cloud/chea3/api/enrich/",
          body = toJSON(body, auto_unbox = TRUE),
          content_type_json(),
          encode = "raw"
        )

        incProgress(0.6, detail = "Processing results...")

        if (http_status(response)$category == "Success") {
          # Parse the response
          results <- content(response, as = "text", encoding = "UTF-8")
          results_list <- fromJSON(results)

          # Store results
          chea3_results(results_list)

          incProgress(1, detail = "Complete")
          showNotification("ChEA3 analysis complete!", type = "message", duration = 3)
        } else {
          showNotification(
            paste("ChEA3 API error:", http_status(response)$message),
            type = "error",
            duration = 5
          )
        }
      }, error = function(e) {
        showNotification(
          paste("Error running ChEA3:", e$message),
          type = "error",
          duration = 5
        )
      })
    })
  })

  # Helper function to render ChEA3 result tables
  render_chea3_table <- function(library_name) {
    req(chea3_results())

    results <- chea3_results()

    # Find the library in results
    lib_data <- NULL
    for (item in results) {
      if (is.data.frame(item) || is.list(item)) {
        if ("Library" %in% names(item) || length(item) > 0) {
          # Check if this is the right library
          if (is.data.frame(item) && nrow(item) > 0) {
            if ("Library" %in% names(item) && any(item$Library == library_name)) {
              lib_data <- item %>% filter(Library == library_name)
              break
            }
          }
        }
      }
    }

    # If not found by Library column, try to find by list name
    if (is.null(lib_data) && library_name %in% names(results)) {
      lib_data <- results[[library_name]]
    }

    # Handle integrated results (Mean Rank)
    if (library_name == "Integrated--meanRank" || library_name == "integrated") {
      for (item in results) {
        if (is.data.frame(item) && "Rank" %in% names(item) && "Score" %in% names(item)) {
          lib_data <- item
          break
        }
      }
    }

    if (!is.null(lib_data) && is.data.frame(lib_data) && nrow(lib_data) > 0) {
      # Select relevant columns if they exist
      cols_to_show <- intersect(names(lib_data), c("Rank", "TF", "Score", "Library", "Overlapping_Genes", "FET p-value", "FDR", "Odds Ratio", "Set_name"))
      if (length(cols_to_show) > 0) {
        lib_data <- lib_data[, cols_to_show, drop = FALSE]
        lib_data$Rank <- as.numeric(lib_data$Rank)

      }

      datatable(
        lib_data,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          searchHighlight = TRUE,
          ordering = TRUE
        ),
        class = 'cell-border stripe hover',
        rownames = FALSE,
        filter = 'top'
      )
    } else {
      datatable(data.frame(Message = "No results available for this library"))
    }
  }

  # Render ChEA3 integrated results
  output$chea3_integrated <- renderDT({
    req(chea3_results())

    results <- chea3_results()

    # The integrated results are typically the first element with meanRank
    integrated_data <- NULL
    for (item in results) {
      if (is.data.frame(item) && "Rank" %in% names(item) && "Score" %in% names(item)) {
        integrated_data <- item
        break
      }
    }
    # Select relevant columns if they exist
    cols_to_show <- intersect(names(integrated_data), c("Rank", "TF", "Score", "Library", "Overlapping_Genes", "FET p-value", "FDR", "Odds Ratio", "Set_name"))
    if (length(cols_to_show) > 0) {
      integrated_data <- integrated_data[, cols_to_show, drop = FALSE]
      integrated_data$Rank <- as.numeric(integrated_data$Rank)
    }

    if (!is.null(integrated_data) && nrow(integrated_data) > 0) {
      datatable(
        integrated_data,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          searchHighlight = TRUE,
          ordering = TRUE
        ),
        class = 'cell-border stripe hover',
        rownames = FALSE,
        filter = 'top'
      )
    } else {
      datatable(data.frame(Message = "No integrated results available"))
    }
  })

  # Render individual library results
  output$chea3_encode <- renderDT({
    req(chea3_results())
    results <- chea3_results()

    lib_data <- NULL
    for (item in results) {
      if (is.data.frame(item) && "Library" %in% names(item)) {
        encode_data <- item %>% filter(grepl("ENCODE", Library, ignore.case = TRUE))
        if (nrow(encode_data) > 0) {
          lib_data <- encode_data
          break
        }
      }
    }

    # Select relevant columns if they exist
    cols_to_show <- intersect(names(lib_data), c("Rank", "TF", "Score", "Library", "Overlapping_Genes", "FET p-value", "FDR", "Odds Ratio", "Set_name"))
    if (length(cols_to_show) > 0) {
      lib_data <- lib_data[, cols_to_show, drop = FALSE]
      lib_data$Rank <- as.numeric(lib_data$Rank)
    }

    if (!is.null(lib_data) && nrow(lib_data) > 0) {
      datatable(lib_data, options = list(pageLength = 15, scrollX = TRUE), class = 'cell-border stripe hover', rownames = FALSE, filter = 'top')
    } else {
      datatable(data.frame(Message = "No ENCODE results available"))
    }
  })

  output$chea3_remap <- renderDT({
    req(chea3_results())
    results <- chea3_results()

    lib_data <- NULL
    for (item in results) {
      if (is.data.frame(item) && "Library" %in% names(item)) {
        remap_data <- item %>% filter(grepl("ReMap", Library, ignore.case = TRUE))
        if (nrow(remap_data) > 0) {
          lib_data <- remap_data
          break
        }
      }
    }

    # Select relevant columns if they exist
    cols_to_show <- intersect(names(lib_data), c("Rank", "TF", "Score", "Library", "Overlapping_Genes", "FET p-value", "FDR", "Odds Ratio", "Set_name"))
    if (length(cols_to_show) > 0) {
      lib_data <- lib_data[, cols_to_show, drop = FALSE]
      lib_data$Rank <- as.numeric(lib_data$Rank)
    }

    if (!is.null(lib_data) && nrow(lib_data) > 0) {
      datatable(lib_data, options = list(pageLength = 15, scrollX = TRUE), class = 'cell-border stripe hover', rownames = FALSE, filter = 'top')
    } else {
      datatable(data.frame(Message = "No ReMap results available"))
    }
  })

  output$chea3_archs4 <- renderDT({
    req(chea3_results())
    results <- chea3_results()

    lib_data <- NULL
    for (item in results) {
      if (is.data.frame(item) && "Library" %in% names(item)) {
        archs4_data <- item %>% filter(grepl("ARCHS4", Library, ignore.case = TRUE))
        if (nrow(archs4_data) > 0) {
          lib_data <- archs4_data
          break
        }
      }
    }

    # Select relevant columns if they exist
    cols_to_show <- intersect(names(lib_data), c("Rank", "TF", "Score", "Library", "Overlapping_Genes", "FET p-value", "FDR", "Odds Ratio", "Set_name"))
    if (length(cols_to_show) > 0) {
      lib_data <- lib_data[, cols_to_show, drop = FALSE]
      lib_data$Rank <- as.numeric(lib_data$Rank)
    }

    if (!is.null(lib_data) && nrow(lib_data) > 0) {
      datatable(lib_data, options = list(pageLength = 15, scrollX = TRUE), class = 'cell-border stripe hover', rownames = FALSE, filter = 'top')
    } else {
      datatable(data.frame(Message = "No ARCHS4 results available"))
    }
  })

  output$chea3_gtex <- renderDT({
    req(chea3_results())
    results <- chea3_results()

    lib_data <- NULL
    for (item in results) {
      if (is.data.frame(item) && "Library" %in% names(item)) {
        gtex_data <- item %>% filter(grepl("GTEx", Library, ignore.case = TRUE))
        if (nrow(gtex_data) > 0) {
          lib_data <- gtex_data
          break
        }
      }
    }

    # Select relevant columns if they exist
    cols_to_show <- intersect(names(lib_data), c("Rank", "TF", "Score", "Library", "Overlapping_Genes", "FET p-value", "FDR", "Odds Ratio", "Set_name"))
    if (length(cols_to_show) > 0) {
      lib_data <- lib_data[, cols_to_show, drop = FALSE]
      lib_data$Rank <- as.numeric(lib_data$Rank)
    }

    if (!is.null(lib_data) && nrow(lib_data) > 0) {
      datatable(lib_data, options = list(pageLength = 15, scrollX = TRUE), class = 'cell-border stripe hover', rownames = FALSE, filter = 'top')
    } else {
      datatable(data.frame(Message = "No GTEx results available"))
    }
  })

  output$chea3_literature <- renderDT({
    req(chea3_results())
    results <- chea3_results()

    lib_data <- NULL
    for (item in results) {
      if (is.data.frame(item) && "Library" %in% names(item)) {
        lit_data <- item %>% filter(grepl("Literature", Library, ignore.case = TRUE))
        if (nrow(lit_data) > 0) {
          lib_data <- lit_data
          break
        }
      }
    }

    # Select relevant columns if they exist
    cols_to_show <- intersect(names(lib_data), c("Rank", "TF", "Score", "Library", "Overlapping_Genes", "FET p-value", "FDR", "Odds Ratio", "Set_name"))
    if (length(cols_to_show) > 0) {
      lib_data <- lib_data[, cols_to_show, drop = FALSE]
      lib_data$Rank <- as.numeric(lib_data$Rank)
    }

    if (!is.null(lib_data) && nrow(lib_data) > 0) {
      datatable(lib_data, options = list(pageLength = 15, scrollX = TRUE), class = 'cell-border stripe hover', rownames = FALSE, filter = 'top')
    } else {
      datatable(data.frame(Message = "No Literature results available"))
    }
  })

  output$chea3_enrichr <- renderDT({
    req(chea3_results())
    results <- chea3_results()

    lib_data <- NULL
    for (item in results) {
      if (is.data.frame(item) && "Library" %in% names(item)) {
        enrichr_data <- item %>% filter(grepl("Enrichr", Library, ignore.case = TRUE))
        if (nrow(enrichr_data) > 0) {
          lib_data <- enrichr_data
          break
        }
      }
    }

    # Select relevant columns if they exist
    cols_to_show <- intersect(names(lib_data), c("Rank", "TF", "Score", "Library", "Overlapping_Genes", "FET p-value", "FDR", "Odds Ratio", "Set_name"))
    if (length(cols_to_show) > 0) {
      lib_data <- lib_data[, cols_to_show, drop = FALSE]
      lib_data$Rank <- as.numeric(lib_data$Rank)
    }

    if (!is.null(lib_data) && nrow(lib_data) > 0) {
      datatable(lib_data, options = list(pageLength = 15, scrollX = TRUE), class = 'cell-border stripe hover', rownames = FALSE, filter = 'top')
    } else {
      datatable(data.frame(Message = "No Enrichr results available"))
    }
  })

  # Generate summary statistics
  output$summary_stats <- renderUI({
    req(filtered_data(), raw())

    # Calculate summary statistics
    filtered <- filtered_data()
    all_data <- raw()
    chea3_data <- chea3_results()

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

    # ChEA3 results
    if (!is.null(chea3_data) && length(chea3_data) > 0) {
      int_chea3 <- NULL
      for (item in chea3_data) {
        if (is.data.frame(item) && "Rank" %in% names(item) && "Score" %in% names(item)) {
          int_chea3 <- item
          break
        }
      } 

      int_chea3$Rank <- as.numeric(int_chea3$Rank)
      top_chea3 <- int_chea3 %>%
        select(Rank, TF, Score, Overlapping_Genes) %>%
        arrange(Rank) %>% 
        head(5)
    }
    else {
      top_chea3 <- NULL
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
          h5(icon("trophy"), "  Top 5 Genes by Effect Size", style = "color: #475569; margin-top: 0;"),
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
          h5(icon("trophy"), "  Top 5 Genes by Mean Effect Score", style = "color: #475569; margin-top: 0;"),
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

      # Row 5: Top TFs table, ChEA3
      if (!is.null(top_chea3) && nrow(top_chea3) > 0) {
        div(
          style = "background: #f8fafc; padding: 15px; border-radius: 6px; margin-bottom: 15px;",
          h5(icon("trophy"), "  Top 5 TFs by ChEA3 Rank", style = "color: #475569; margin-top: 0;"),
          h5("    (based on genes and pathway selected in ChEA3 TF Analysis tab)", style = "color: #94a3b8; font-style: italic; margin-top: 0; font-size: 12px;"),
          tags$table(
            class = "table table-striped table-sm",
            style = "font-size: 13px; margin-bottom: 0;",
            tags$thead(
              tags$tr(
                tags$th("Rank"),
                tags$th("TF"),
                tags$th("Score"),
                tags$th("Overlapping Genes")
              )
            ),
            tags$tbody(
              lapply(1:nrow(top_chea3), function(i) {
                tags$tr(
                  tags$td(strong(top_chea3$Rank[i])),
                  tags$td(strong(top_chea3$TF[i])),
                  tags$td(strong(top_chea3$Score[i])),
                  tags$td(strong(top_chea3$Overlapping_Genes[i]))
                )
              })
            )
          )
        )
      } 
      else {
         p("Please run ChEA3 for top TFs", style = "color: #94a3b8; font-style: italic;")
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

  # Download ChEA3 results
  output$download_chea3 <- downloadHandler(
    filename = function() {
      paste0("chea3_results_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      req(chea3_results())

      wb <- createWorkbook()
      results <- chea3_results()

      # Add each result set as a sheet
      sheet_num <- 1
      for (i in seq_along(results)) {
        item <- results[[i]]
        if (is.data.frame(item) && nrow(item) > 0) {
          sheet_name <- names(results)[i]
          if (is.null(sheet_name) || sheet_name == "") {
            sheet_name <- paste0("Results_", sheet_num)
          }
          # Truncate sheet name to 31 characters (Excel limit)
          sheet_name <- substr(sheet_name, 1, 31)
          addWorksheet(wb, sheet_name)
          writeData(wb, sheet_name, item)
          sheet_num <- sheet_num + 1
        }
      }

      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

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
