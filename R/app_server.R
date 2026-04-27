#' @import shiny
#' @import shinythemes
#' @import shinyjs
#' @import viridisLite
#' @import DT
#' @import shinyWidgets
#' @import shinyBS
#' @import shinycssloaders
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @import data.table
#' @import stringr
#' @import vroom
#' @import plotly
app_server <- function(input, output, session) {

  session$onFlushed(function() {
    shinyjs::disable("run_proc")
    shinyjs::disable("dl_volcano")
    shinyjs::disable("dl_matrix")
  }, once = TRUE)

  observe({
    shinyjs::toggleState("run_proc", condition = !is.null(input$file_data))
  })

  observe({
    if (procReady()) {
      shinyjs::enable("dl_volcano")
      shinyjs::enable("dl_matrix")
    } else {
      shinyjs::disable("dl_volcano")
      shinyjs::disable("dl_matrix")
    }
  })

  observeEvent(input$file_data, {
    rv$raw <- NULL
    rv$mat <- NULL
    rv$fmap <- NULL
    rv$labels <- NULL
    rv$df_used <- NULL
    rv$volcano <- NULL
  }, ignoreInit = TRUE)

  rv <- reactiveValues(
    raw = NULL,
    mat = NULL,
    fmap = NULL,
    labels = NULL,
    df_used = NULL,
    volcano = NULL
  )

  procReady <- reactive({
    !is.null(rv$volcano) && nrow(rv$volcano) > 0
  })

  # ---- Load raw data (Robust Switch) ----
  raw_df <- reactive({
    req(input$file_data)
    ext <- tools::file_ext(input$file_data$name)
    validate(need(tolower(ext) == "csv", "Please upload a .csv file"))

    tool <- input$software_tool %||% "mzmine"

    if (tool == "msdial") {
      # Use custom robust reader for MS-DIAL
      df <- read_msdial_robust(input$file_data$datapath)
      clean_mzmine_export(df) # cleanup empty trailing cols
    } else {
      # Use fast vroom reader for standard CSVs (xcms/mzmine/default)
      clean_mzmine_export(vroom::vroom(input$file_data$datapath, delim = ","))
    }
  }) %>% bindCache(input$file_data$name, input$software_tool)

  # ---- Column Pickers (Smart Defaults) ----
  output$col_pickers <- renderUI({
    req(raw_df())
    cols <- names(raw_df())
    tool <- input$software_tool %||% "mzmine"

    rid_cand <- switch(tool,
                       xcms   = c("...1", "X...1", "feature_id", "feature", "id", "row id"),
                       msdial = c("alignment id", "alignmentid", "spot id"),
                       default = c("Feature", "row id", "id", "feature_id"),
                       c("row id", "id", "feature_id")
    )
    mz_cand <- switch(tool,
                      xcms   = c("mzmed", "mz", "m/z", "mzmin", "mzmax"),
                      msdial = c("average mz", "averagemz", "mz"),
                      default = c("mz", "m/z", "mass", "average mz"),
                      c("row m/z", "row mz", "mz")
    )
    rt_cand <- switch(tool,
                      xcms   = c("rtmed", "rt", "rtmin"),
                      msdial = c("average rt(min)", "average rt", "averagertmin", "rt"),
                      default = c("rt", "retention time", "time"),
                      c("row retention time", "row rt", "rt")
    )

    def_rid <- guess_col(cols, rid_cand)
    def_mz  <- guess_col(cols, mz_cand)
    def_rt  <- guess_col(cols, rt_cand)

    choices_rid <- c("<Auto-generate>", cols)
    sel_rid <- if (is.null(def_rid)) "<Auto-generate>" else def_rid

    tagList(
      selectInput("row_id_col", "Row ID column:", choices = choices_rid, selected = sel_rid),
      selectInput("mz_col",     "m/z column:",   choices = cols, selected = def_mz  %||% cols[1]),
      selectInput("rt_col",     "rt column:",    choices = cols, selected = def_rt  %||% cols[1])
    )
  })

  output$raw_header <- renderUI({
    req(raw_df())
    h3(sprintf("Raw dataset: %d rows × %d columns", nrow(raw_df()), ncol(raw_df())))
  })

  output$raw_preview <- renderDT({
    req(raw_df())
    datatable(head(raw_df(), 20), options = list(scrollX = TRUE, pageLength = 8))
  })

  # ---- Build matrix + fmap ----
  built <- reactive({
    req(raw_df(), input$row_id_col, input$mz_col, input$rt_col)

    df <- raw_df()
    rid_col <- input$row_id_col

    if (rid_col == "<Auto-generate>") {
      df$FeatureID_Auto <- as.numeric(seq_len(nrow(df)))
      rid_col <- "FeatureID_Auto"
    }

    kws <- input$sample_keywords
    if (is.null(kws) || length(kws) == 0) kws <- c(".")

    parse_feature_table_to_matrix(
      df,
      row_id_col = rid_col,
      mz_col     = input$mz_col,
      rt_col     = input$rt_col,
      sample_keywords = kws
    )
  })

  sample_names <- reactive({
    req(built())
    rownames(built()$mat)
  })

  # ---- Labels ----
  labels_vec <- reactive({
    req(sample_names())

    if (input$label_source == "csv") {
      req(input$file_labels)
      v <- read_onecol_csv(input$file_labels$datapath)
      v <- trimws(v)
      validate(need(length(v) == length(sample_names()),
                    sprintf("Labels count (%d) must match #samples (%d).",
                            length(v), length(sample_names()))))
      v
    } else {
      sn <- sample_names()
      sn2 <- if (isTRUE(input$clean_sample_names)) clean_sample_names(sn) else sn
      sep <- input$token_sep %||% "_"
      idx <- as.integer(input$token_index %||% 2)

      parts <- strsplit(sn2, sep, fixed = TRUE)
      ok <- vapply(parts, function(z) length(z) >= idx, logical(1))
      validate(need(all(ok), sprintf("Token %d missing in some sample names. Adjust separator/index.", idx)))
      labs <- vapply(parts, function(z) z[[idx]], character(1))
      validate(need(all(nzchar(labs)), "Parsed empty labels — adjust separator/index."))
      labs
    }
  })

  output$labels_header <- renderUI({
    req(sample_names(), labels_vec())
    h3(sprintf("Labels ready: %d samples", length(labels_vec())))
  })

  output$labels_table <- renderDT({
    req(sample_names(), labels_vec())
    datatable(
      tibble(Sample = sample_names(), Label = labels_vec()),
      options = list(pageLength = 8, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$ref_group_picker <- renderUI({
  req(labels_vec())
  levs <- sort(unique(labels_vec()))
  selectInput("ref_group", "Reference Group:", choices = levs)
  })

  # ---- SIRIUS pickers ----
  sirius_df <- reactive({
    req(input$use_sirius)
    req(input$file_sirius)
    ext <- tools::file_ext(input$file_sirius$name)
    validate(need(tolower(ext) == "csv", "SIRIUS file must be .csv"))
    vroom::vroom(input$file_sirius$datapath, delim = ",") %>% as.data.frame(check.names = FALSE)
  })

  output$sirius_pickers <- renderUI({
    req(sirius_df())
    cols <- names(sirius_df())
    tagList(
      selectInput("sirius_idcol", "SIRIUS Feature ID column:",
                  choices = cols, selected = if ("mappingFeatureId" %in% cols) "mappingFeatureId" else cols[1]),
      selectInput("sirius_npcol", "NPC column:",
                  choices = cols, selected = if ("NPC#class" %in% cols) "NPC#class" else cols[1]),
      selectInput("sirius_cfcol", "ClassyFire column:",
                  choices = cols, selected = if ("ClassyFire#class" %in% cols) "ClassyFire#class" else cols[1])
    )
  })

  # ---- Process button ----
  observeEvent(input$run_proc, {
    req(built(), labels_vec())
    labs_pre <- labels_vec()
    if (stop_if_one_group(labs_pre)) return(NULL)
    withProgress(message = "Processing...", value = 0, {
      incProgress(0.15, detail = "Building matrix")
      X <- built()$mat
      fmap <- built()$fmap

      incProgress(0.20, detail = "Adding labels")
      labs <- labels_vec()
      validate(need(length(labs) == nrow(X),
                    sprintf("Labels length (%d) must match #samples (%d).",
                            length(labs), nrow(X))))

      df_raw <- as.data.frame(X, check.names = FALSE, stringsAsFactors = FALSE)
      df_used <- cbind(Label = labs, df_raw)
      df_used$Label <- as.factor(df_used$Label)

      incProgress(0.25, detail = "Imputation (if enabled)")
      if (identical(input$do_mvi, "yes")) {
          X0 <- df_used[, -1, drop = FALSE]
          Xm <- impute_lod_random(
            X0,
            noise_mode     = input$noise_mode %||% "quantile",
            noise_quantile = input$noise_quantile %||% 0.25,
            noise_manual   = input$noise_manual %||% 50,
            sd_val         = input$noise_sd %||% 30,
            seed           = 1234
          )
        df_used <- as.data.frame(cbind(Label = df_used$Label, as.data.frame(Xm, check.names = FALSE)),
                                 check.names = FALSE, stringsAsFactors = FALSE)
        df_used$Label <- as.factor(df_used$Label)
      }

      incProgress(0.30, detail = "Running...")
      volc <- compute_stats_long(
      df_used,
      test = input$test_type %||% "Student",
      adj  = input$p_adjust %||% "BH",
      paired = isTRUE(input$paired),
      eqvar  = isTRUE(input$eqvar),
      pseudocount = 1.1,
      log2_test = isTRUE(input$log2_test),
      ref_group = input$ref_group
      )

      incProgress(0.05, detail = "Joining mz/rt/id")
      volc <- volc %>% left_join(fmap, by = "Feature")

      # default annotation columns
      volc$`NPC#class` <- "Not provided"
      volc$`ClassyFire#class` <- "Not provided"

      incProgress(0.05, detail = "Joining SIRIUS (optional)")
      if (isTRUE(input$use_sirius) && !is.null(input$file_sirius)) {
        s <- sirius_df()
        req(input$sirius_idcol, input$sirius_npcol, input$sirius_cfcol)

        ss <- s %>%
          transmute(
            id = as.character(.data[[input$sirius_idcol]]),
            `NPC#class` = as.character(.data[[input$sirius_npcol]]),
            `ClassyFire#class` = as.character(.data[[input$sirius_cfcol]])
          )

        volc <- volc %>%
          mutate(id = as.character(id)) %>%
          left_join(ss, by = "id", suffix = c("", ".sirius")) %>%
          mutate(
            `NPC#class` = dplyr::coalesce(.data[["NPC#class.sirius"]], .data[["NPC#class"]]),
            `ClassyFire#class` = dplyr::coalesce(.data[["ClassyFire#class.sirius"]], .data[["ClassyFire#class"]])
          ) %>%
          select(-any_of(c("NPC#class.sirius", "ClassyFire#class.sirius")))
      }

      volc <- volc %>%
        mutate(
          mz = round(mz, 6),
          RT = round(RT, 4),
          Mean = signif(Mean, 4),
          FC = round(FC, 4),
          `Adj.p-value` = as.numeric(`Adj.p-value`),
          `Adj.p-value.log` = round(`Adj.p-value.log`, 5)
        )

      rv$raw <- built()$raw
      rv$mat <- X
      rv$fmap <- fmap
      rv$labels <- labs
      rv$df_used <- df_used
      rv$volcano <- volc
    })

    showNotification("Processing finished. Switching to Volcano explorer...", type = "message", duration = 4)
    updateTabsetPanel(session, "tabs", selected = "volcano")
  }, ignoreInit = TRUE)

  output$proc_summary <- renderUI({
    if (!procReady()) {
      div(class = "highlight", "Not processed yet. Upload -> Set settings -> click 'Run preprocessing'.")
    } else {
      div(
        style="background:#ffffffcc; padding:10px; border-radius:10px; font-size:20px; border:1px solid #ddd;",
        h4("Processing summary"),
        tags$ul(
          tags$li(sprintf("Samples: %d", nrow(rv$df_used))),
          tags$li(sprintf("Features: %d", ncol(rv$df_used) - 1)),
          tags$li(sprintf("Comparisons: %d", length(unique(rv$volcano$Groups)))),
          tags$li(sprintf("Rows in volcano table: %d", nrow(rv$volcano)))
        ),
        div(class="small-note", "")
      )
    }
  })

  # ---------------- Volcano tab UI ----------------

      observeEvent(input$reset_filters, {
      req(procReady())

      # switches
      updateMaterialSwitch(session, "sig_only", value = FALSE)
      updateMaterialSwitch(session, "use_fdr_filter", value = FALSE)
      updateMaterialSwitch(session, "use_fc_filter",  value = FALSE)

      # pickers / radios
      updatePickerInput(session, "sel_feat", selected = character(0))
      updateRadioButtons(session, "color_by", selected = "Groups")
      updateRadioButtons(session, "present_as", selected = "Boxplot")
      updateRadioButtons(session, "fc_dir", selected = "both")

      # sliders
      dd <- rv$volcano
      mzr <- finite_range(dd$mz)
      rtr <- finite_range(dd$RT)
      mr  <- finite_range(log10(dd$Mean + 1.1))

      if (!is.null(mzr)) updateSliderInput(session, "mz_range", value = c(mzr[1], mzr[2]))
      if (!is.null(rtr)) updateSliderInput(session, "rt_range", value = c(rtr[1], rtr[2]))
      if (!is.null(mr))  updateSliderInput(session, "intensity_range", value = c(round(mr[1], 1), round(mr[2], 1)))

      updateSliderInput(session, "fdr_cutoff", value = 1.30103)
      updateSliderInput(session, "fc_thr", value = 1)
    })

  output$volcano_sidebar <- renderUI({
    if (!procReady()) {
      return(div(class="highlight",
                 "No processed dataset yet. Go to '1) Load & Process' and click 'Run preprocessing'."))
    }

    tagList(
      actionButton("reset_filters", "Reset filters", class = "btn btn-warning"),
      br(),
      br(),
      materialSwitch("sig_only", "Significant only", value = FALSE, status = "info"),
        pickerInput(
        inputId = "sel_feat",
        label   = "Select/deselect features (optional):",
        choices = sort(unique(rv$volcano$Feature)),
        options = list(`actions-box` = TRUE, `live-search` = TRUE),
        multiple = TRUE
      ),
      br(),
      radioButtons(
        "color_by",
        "Color points by:",
        choices = c("Groups" = "Groups", "Mean" = "Mean"),
        selected = "Groups",
        inline = TRUE
      ),

      radioButtons(
        "present_as",
        "Click plot shows:",
        choices = c("Boxplot" = "Boxplot", "Scatterplot" = "Scatterplot"),
        selected = "Boxplot",
        inline = TRUE
      ),

      uiOutput("volcano_sliders")
    )
  })

  output$volcano_main <- renderUI({
    if (!procReady()) return(NULL)
    tagList(
      withSpinner(plotlyOutput("volcano_plot", height = "520px"), type = 8, color = "#66CDAA"),
      div(style = "height:8px;"),
      plotlyOutput("feature_plot", height = "260px")
    )
  })

  output$volcano_sliders <- renderUI({
    req(procReady())
    dd <- rv$volcano

    mzr <- finite_range(dd$mz)
    rtr <- finite_range(dd$RT)
    mr  <- finite_range(log10(dd$Mean+ 1.1))
    fcmax <- max(abs(dd$FC), na.rm = TRUE)
    yMax <- max(dd$`Adj.p-value.log`, na.rm = TRUE)

    validate(need(!is.null(mzr) && !is.null(rtr) && !is.null(mr),
                  "No finite mz/RT/Mean values available for sliders."))

    tagList(
      sliderInput("mz_range", "m/z:",
                  min = round(min(dd$mz, na.rm = TRUE), 3),
                  max = round(max(dd$mz, na.rm = TRUE), 3),
                  value = round(range(dd$mz, na.rm = TRUE), 3),
                  step = 0.1),
      sliderInput("rt_range", "RT:",
                  min = round(min(dd$RT, na.rm = TRUE), 1),
                  max = round(max(dd$RT, na.rm = TRUE), 1),
                  value = round(range(dd$RT, na.rm = TRUE),1),
                  step = 0.1),
      sliderInput("intensity_range", "Mean log10(Intensity):",
                  min = round(log10(min(dd$Mean+1.1, na.rm = TRUE)),1),
                  max = round(log10(max(dd$Mean+1.1, na.rm = TRUE)),1),
                  value = round(log10(range(dd$Mean+1.1, na.rm = TRUE)),1),
                  step = 0.1),
          tags$hr(),
       materialSwitch("use_fdr_filter", "Filter by Adj.p-value", value = FALSE, status = "success"),
    conditionalPanel(
      condition = "input.use_fdr_filter == true",
      sliderInput("fdr_cutoff", "-log10(Adj.p-value):",
                  min = 0, max = yMax, value = 1.30103, step = 0.01),
      uiOutput("fdr_equiv_text")
    ),

    materialSwitch("use_fc_filter", "Filter by Fold-Change", value = FALSE, status = "success"),
    conditionalPanel(
      condition = "input.use_fc_filter == true",
      sliderInput("fc_thr", "FC threshold (|log2FC| ≥):",
                  min = 0, max = round(fcmax, 1), value = 1, step = 0.1),
      radioButtons(
        "fc_dir",
        "Direction:",
        choices = c("Both sides" = "both", "Up only" = "up", "Down only" = "down"),
        selected = "both",
        inline = TRUE
      )
    ))
  })

  # ---- Filtered volcano data
  filtered_volcano <- reactive({
    req(procReady(), input$mz_range, input$rt_range, input$intensity_range)

    dd <- rv$volcano

    if (!is.null(input$sel_feat) && length(input$sel_feat) > 0) {
      return(dd %>% dplyr::filter(Feature %in% input$sel_feat))
    }
    if (isTRUE(input$sig_only)) {
      dd <- dd %>% dplyr::filter(Significant)
    }

    dd <- dd %>%
      dplyr::filter(
        is.finite(mz), is.finite(RT), is.finite(Mean),
        mz >= input$mz_range[1], mz <= input$mz_range[2],
        RT >= input$rt_range[1], RT <= input$rt_range[2],
        log10(Mean + 1.1) >= input$intensity_range[1],
        log10(Mean + 1.1) <= input$intensity_range[2]
      )

    if (isTRUE(input$use_fdr_filter)) {
      thr <- input$fdr_cutoff %||% 1.30103
      dd <- dd %>% dplyr::filter(is.finite(`Adj.p-value.log`), `Adj.p-value.log` >= thr)
    }

    output$fdr_equiv_text <- renderUI({
      req(input$fdr_cutoff)
      tags$small(sprintf("Equivalent Adj.p-value ≤ %.3g", 10^(-input$fdr_cutoff)))
    })

    if (isTRUE(input$use_fc_filter)) {
      thr <- input$fc_thr
      if (input$fc_dir == "both") {
        dd <- dd %>% dplyr::filter(abs(FC) >= thr)
      } else if (input$fc_dir == "up") {
        dd <- dd %>% dplyr::filter(FC >= thr)
      } else {
        dd <- dd %>% dplyr::filter(FC <= -thr)
      }
    }

    dd
  })

  # ---- Volcano plot
  output$volcano_plot <- renderPlotly({
    req(filtered_volcano())
    dd <- filtered_volcano()
    validate(need(nrow(dd) > 0, "No points left after filtering."))

    dd$key <- paste(dd$Groups, dd$Feature, sep = "__")

    hover_txt <- paste0(
      "Groups: ", dd$Groups,
      "<br>FC: ", dd$FC,
      "<br>FDR: ", format(dd$`Adj.p-value`, digits = 3, scientific = TRUE),
      "<br>Test scale: ", dd$TestScale,
      "<br>Feature: ", dd$Feature,
      "<br>ID: ", dd$id,
      "<br>m/z: ", dd$mz,
      "<br>RT: ", dd$RT,
      "<br>Mean Intensity: ", format(dd$Mean, big.mark = ",", scientific = FALSE),
      "<br>NPC: ", dd$`NPC#class`,
      "<br>ClassyFire: ", dd$`ClassyFire#class`
    )

    fc_line <- if (isTRUE(input$use_fc_filter)) input$fc_thr else 1
    ythr <- if (isTRUE(input$use_fdr_filter)) (input$fdr_cutoff %||% 1.30103) else 1.30103

    shapes <- list()

    if (isTRUE(input$use_fc_filter)) {
      if (input$fc_dir %in% c("both", "down")) {
        shapes <- c(shapes, list(list(type="line", x0=-fc_line, x1=-fc_line, xref="x",
                                      y0=0, y1=1, yref="paper", line=list(dash="dot"))))
      }
      if (input$fc_dir %in% c("both", "up")) {
        shapes <- c(shapes, list(list(type="line", x0= fc_line, x1= fc_line, xref="x",
                                      y0=0, y1=1, yref="paper", line=list(dash="dot"))))
      }
    } else {
      shapes <- c(shapes, list(
        list(type="line", x0=-1, x1=-1, xref="x", y0=0, y1=1, yref="paper", line=list(dash="dot")),
        list(type="line", x0= 1, x1= 1, xref="x", y0=0, y1=1, yref="paper", line=list(dash="dot"))
      ))
    }

    shapes <- c(shapes, list(
      list(type="line", x0=0, x1=1, xref="paper", y0=ythr, y1=ythr, yref="y", line=list(dash="dot"))
    ))

    if (input$color_by == "Groups") {
      plot_ly(
        data = dd,
        x = ~FC, y = ~`Adj.p-value.log`,
        colors = "Set1",
        color = ~Groups,
        type = "scatter", mode = "markers",
        text = hover_txt, hoverinfo = "text",
        key = ~key,
        marker = list(size = 12, opacity = 0.85, line = list(color = "black", width = 1)),
        source = "volcano"
      ) %>%
        layout(
          shapes = shapes,
          xaxis = list(title = "log2(FC)"),
          yaxis = list(title = "-log10(FDR)"),
          legend = list(title = list(text = "Comparison"))
        ) %>% event_register("plotly_click")
    } else {
      plot_ly(
        data = dd,
        x = ~FC, y = ~`Adj.p-value.log`,
        color = ~log10(Mean + 1.1),
        symbol = ~Groups,
        colors = viridis(12),
        type = "scatter", mode = "markers",
        text = hover_txt, hoverinfo = "text",
        key = ~key,
        marker = list(size = 12, opacity = 0.9, line = list(color = "black", width = 1)),
        source = "volcano"
      ) %>%
        layout(
          shapes = shapes,
          xaxis = list(title = "log2(FC)"),
          yaxis = list(title = "-log10(FDR)"),
          legend = list(title = list(text = "Comparison"))
        ) %>% event_register("plotly_click")
    }
  })

  # ---- Feature plot on click
  output$feature_plot <- renderPlotly({
    req(procReady(), rv$df_used, rv$mat)

    click <- event_data("plotly_click", source = "volcano")
    if (is.null(click) || is.null(click$key)) return(NULL)

    key <- click$key[[1]]
    parts <- str_split_fixed(key, "__", 2)
    comp  <- parts[1, 1]
    feat  <- parts[1, 2]

    df_used <- rv$df_used
    validate(need(feat %in% colnames(df_used), "Clicked feature not found in matrix."))

    row <- rv$volcano %>% filter(Groups == comp, Feature == feat) %>% slice(1)
    s_names <- rownames(rv$mat)

    yy <- as.numeric(df_used[[feat]])
    xx <- as.character(df_used$Label)

    if (input$present_as == "Boxplot") {
      plot_ly(
        x = xx, y = yy,
        colors = "Dark2",
        type = "box",
        color = xx,
        boxpoints = FALSE,
        marker = list(opacity = 0.8)
      ) %>% hide_legend() %>%
        layout(
          title = list(text = paste0(feat, "<br><span style='font-size:12px;'>")),
          xaxis = list(title = "Group"),
          yaxis = list(title = "Intensity")
        )
    } else {
      plot_ly(
        x = xx, y = yy,
        text = s_names,
        colors = "Dark2",
        type = "scatter", mode = "markers",
        color = xx,
        hovertemplate = paste(
          "<b>Sample:</b> %{text}<br>",
          "<b>Group:</b> %{x}<br>",
          "<b>Intensity:</b> %{y}",
          "<extra></extra>"
        ),
        marker = list(size = 20, opacity = 0.85, symbol = "diamond", line = list(color = "black", width = 2))
      ) %>% hide_legend() %>%
        layout(
          title = list(text = paste0(feat, "<br><span style='font-size:12px;'>")),
          xaxis = list(title = "Group"),
          yaxis = list(title = "Intensity")
        )
    }
  })

  # ---- Downloads
  output$dl_volcano <- downloadHandler(
    filename = function() "volcano_table.csv",
    content = function(file) {
      req(rv$volcano)
      data.table::fwrite(rv$volcano, file)
    }
  )

output$dl_matrix <- downloadHandler(
  filename = function() "processed_table.csv",
  content = function(file) {
    req(rv$df_used, rv$mat)

    out <- as.data.frame(rv$df_used, check.names = FALSE, stringsAsFactors = FALSE)
    out <- cbind(Sample = rownames(rv$mat), out)

    data.table::fwrite(out, file)
  }
)
}
