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
#' @import limma
#' @import RColorBrewer
#' @import zip
app_server <- function(input, output, session) {

  output$label_upload_warning <- renderUI({
  src <- input$label_source %||% "token"

  need_file <- FALSE

  if (identical(src, "csv") && is.null(input$file_labels)) {
    need_file <- TRUE
  }

  if (identical(src, "metadata") && is.null(input$file_metadata_labels)) {
    need_file <- TRUE
  }

  if (!need_file) return(NULL)

  div(
    class = "alert alert-warning",
    style = "
      margin-top: 10px;
      margin-bottom: 12px;
      font-size: 16px;
      font-weight: 700;
      border: 2px solid #f0ad4e;
      border-radius: 8px;
    ",
    "Please upload a CSV file first!"
  )
})

  session$onFlushed(function() {
    shinyjs::disable("run_proc")
    shinyjs::disable("dl_volcano")
    shinyjs::disable("dl_matrix")
    shinyjs::disable("dl_autoplotter_zip")
  }, once = TRUE)

  observe({
    shinyjs::toggleState("run_proc", condition = !is.null(input$file_data))
  })

 observe({
  ids <- c(
    "dl_volcano",
    "dl_matrix",
    "dl_autoplotter_zip"
  )

  if (procReady()) {
    lapply(ids, shinyjs::enable)
  } else {
    lapply(ids, shinyjs::disable)
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

  observeEvent(
  list(
    input$label_source,
    input$file_labels,
    input$file_metadata_labels,
    input$metadata_sample_col,
    input$metadata_label_col,
    input$metadata_clean_sample_names,
    input$metadata_remove_suffixes,
    input$token_sep,
    input$token_index
  ),
  {
    rv$labels <- NULL
    rv$df_used <- NULL
    rv$volcano <- NULL
  },
  ignoreInit = TRUE
)

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

  dataset_name <- reactive({
  nm <- input$file_data$name %||% "dataset.csv"
  tools::file_path_sans_ext(basename(nm))
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

  manual_labels <- reactiveVal(NULL)

auto_label_table <- reactive({
  req(sample_names())

  make_label_table(
    sample_names(),
    labels_from_sample_names_or_raw(
      sample_names(),
      token_sep = input$token_sep %||% "_",
      token_index = input$token_index %||% 2,
      clean_names = FALSE
    )
  )
})

observeEvent(sample_names(), {
  req(auto_label_table())
  manual_labels(auto_label_table())
}, ignoreInit = FALSE)

observeEvent(input$fill_manual_labels, {
  req(auto_label_table())

  manual_labels(auto_label_table())

  rv$labels <- NULL
  rv$df_used <- NULL
  rv$volcano <- NULL

  showNotification(
    "Editable label table was filled from current token labels.",
    type = "message",
    duration = 3
  )
}, ignoreInit = TRUE)

observeEvent(input$labels_table_cell_edit, {
  req(input$label_source == "manual")

  info <- input$labels_table_cell_edit

  tbl <- manual_labels()
  req(tbl)

  row_i <- as.integer(info$row)

  if (!is.finite(row_i) || row_i < 1 || row_i > nrow(tbl)) {
    showNotification("Edited row is outside label table.", type = "error", duration = 3)
    return(NULL)
  }

  tbl$Label[row_i] <- trimws(as.character(info$value))

  manual_labels(tbl)

  rv$labels <- NULL
  rv$df_used <- NULL
  rv$volcano <- NULL

  showNotification(
    paste0("Label updated: ", tbl$Sample[row_i], " -> ", tbl$Label[row_i]),
    type = "message",
    duration = 2
  )
}, ignoreInit = TRUE)

metadata_raw <- reactive({
  req(input$file_metadata_labels)
  read_metadata_csv(input$file_metadata_labels, context = "metadata labels")
})


output$metadata_sample_col_ui <- renderUI({
  req(metadata_raw())

  cols <- names(metadata_raw())

  selectInput(
    "metadata_sample_col",
    "Metadata sample-name column:",
    choices = cols,
    selected = guess_metadata_sample_col(cols)
  )
})


output$metadata_label_col_ui <- renderUI({
  req(metadata_raw())

  cols <- names(metadata_raw())
  sample_col <- input$metadata_sample_col %||% guess_metadata_sample_col(cols)
  choices <- setdiff(cols, sample_col)

  validate(
    need(length(choices) > 0, "Metadata file has no column available for labels.")
  )

  selectizeInput(
    "metadata_label_col",
    "Metadata column to use as Label:",
    choices = choices,
    selected = guess_metadata_label_col(cols, sample_col),
    multiple = FALSE
  )
})


metadata_labels <- reactive({
  req(
    input$file_metadata_labels,
    input$metadata_sample_col,
    input$metadata_label_col,
    sample_names()
  )

  metadata_labels_by_sample(
    upload = input$file_metadata_labels,
    sample_names = sample_names(),
    sample_col = input$metadata_sample_col,
    label_col = input$metadata_label_col,
    clean_enabled = isTRUE(input$metadata_clean_sample_names),
    remove_suffixes = input$metadata_remove_suffixes %||% character(0),
    context = "metadata labels"
  )
})

  # ---- Labels ----
  labels_vec <- reactive({
  req(sample_names())

  src <- input$label_source %||% "token"

  if (identical(src, "csv")) {

  req(input$file_labels)

  v <- read_onecol_csv(input$file_labels$datapath)
  v <- trimws(v)

  validate(
    need(
      length(v) == length(sample_names()),
      sprintf(
        "Labels count (%d) must match #samples (%d).",
        length(v), length(sample_names())
      )
    )
  )

  v

} else if (identical(src, "metadata")) {

  metadata_labels()

} else if (identical(src, "manual")) {

    tbl <- manual_labels()
    req(tbl)

    validate(
      need(
        nrow(tbl) == length(sample_names()),
        "Manual label table must match the number of samples."
      ),
      need(
        identical(as.character(tbl$Sample), as.character(sample_names())),
        "Manual label table does not match current sample names. Click 'Fill editable table from current token labels'."
      ),
      need(
        !any(is.na(tbl$Label) | trimws(tbl$Label) == ""),
        "All samples must have labels."
      )
    )

    trimws(as.character(tbl$Label))

  } else {

    sn <- sample_names()
    sn2 <- sn
    sep <- input$token_sep %||% "_"
    idx <- as.integer(input$token_index %||% 2)

    parts <- strsplit(sn2, sep, fixed = TRUE)
    ok <- vapply(parts, function(z) length(z) >= idx, logical(1))

    validate(
      need(
        all(ok),
        sprintf("Token %d missing in some sample names. Adjust separator/index.", idx)
      )
    )

    labs <- vapply(parts, function(z) z[[idx]], character(1))

    validate(
      need(all(nzchar(labs)), "Parsed empty labels — adjust separator/index.")
    )

    labs
  }
})

  output$labels_header <- renderUI({
    req(sample_names(), labels_vec())
    h3(sprintf("Labels ready: %d samples", length(labels_vec())))
  })

  output$labels_table <- renderDT({
  req(sample_names())

  src <- input$label_source %||% "token"

  tbl <- if (identical(src, "manual")) {
    manual_labels()
  } else {
    req(labels_vec())
    make_label_table(sample_names(), labels_vec())
  }

  req(tbl)

  datatable(
    tbl,
    editable = if (identical(src, "manual")) {
      list(
        target = "cell",
        disable = list(columns = c(0)) # lock Sample column
      )
    } else {
      FALSE
    },
    options = list(
      pageLength = 8,
      scrollX = TRUE,
      ordering = FALSE,
      searching = FALSE
    ),
    rownames = FALSE
  )
}, server = FALSE)

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

  # ---- SIRIUS & GNPS stats tab ----

guess_col_ci <- function(cols, candidates, default = NULL) {
  if (is.null(cols) || !length(cols)) return(default)

  cols_norm <- gsub("[^a-z0-9]", "", tolower(cols))
  cand_norm <- gsub("[^a-z0-9]", "", tolower(candidates))

  for (cand in cand_norm) {
    hit <- which(cols_norm == cand)
    if (length(hit)) return(cols[hit[1]])
  }

  for (cand in cand_norm) {
    hit <- which(grepl(cand, cols_norm, fixed = TRUE))
    if (length(hit)) return(cols[hit[1]])
  }

  if (!is.null(default)) default else cols[1]
}

clean_stats_value <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x) | x %in% c("NA", "NaN", "null", "Not provided")] <- NA_character_
  x
}

gnps_pairs_df <- reactive({
  req(input$file_gnps_pairs)

  ext <- tolower(tools::file_ext(input$file_gnps_pairs$name))

  validate(
    need(
      ext %in% c("tsv", "txt", "csv"),
      "Upload GNPS pairs as .tsv, .txt, or .csv."
    )
  )

  delim <- if (ext == "csv") "," else "\t"

  as.data.frame(
    vroom::vroom(
      input$file_gnps_pairs$datapath,
      delim = delim,
      col_names = TRUE,
      show_col_types = FALSE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})

output$sirius_gnps_sidebar <- renderUI({
  if (!isTRUE(input$use_sirius) || is.null(input$file_sirius)) {
    return(
      div(
        class = "highlight",
        "Switch on 'Join SIRIUS summary' in the Load & Process tab and upload the SIRIUS .csv file."
      )
    )
  }

  s <- sirius_df()
  s_cols <- names(s)

  default_id <- if (!is.null(input$sirius_idcol) && input$sirius_idcol %in% s_cols) {
    input$sirius_idcol
  } else {
    guess_col_ci(
      s_cols,
      c("mappingFeatureId", "id", "featureId", "row ID"),
      s_cols[1]
    )
  }

  default_ann <- guess_col_ci(
    s_cols,
    c("NPC#class", "ClassyFire#class", "NPC class", "ClassyFire class", "class", "superclass"),
    s_cols[1]
  )

  tagList(
    h3(class = "highlight", "SIRIUS annotation frequency"),

    selectInput(
      "stats_sirius_id_col",
      "SIRIUS ID column:",
      choices = s_cols,
      selected = default_id
    ),

    selectInput(
      "stats_sirius_col",
      "SIRIUS column for frequency statistics:",
      choices = s_cols,
      selected = default_ann
    ),

    uiOutput("stats_class_picker"),

    tags$hr(),

    h3(class = "highlight", "Optional GNPS network pairs"),

    fileInput(
      "file_gnps_pairs",
      "Upload GNPS Pairs List file (.tsv/.txt/.csv)",
      accept = c(".tsv", ".txt", ".csv")
    ),

    uiOutput("gnps_pickers_ui"),

    div(
      class = "small-note",
      "GNPS pairs are converted from ClusterID1/ClusterID2 into one ClusterID column. By default, no peak-table pruning is applied. If a peak-table ID column is selected, SIRIUS IDs and GNPS ClusterIDs are restricted to IDs present in the uploaded peak table."
    )
  )
})

output$gnps_pickers_ui <- renderUI({
  req(gnps_pairs_df())

  g_cols <- names(gnps_pairs_df())

  raw_cols <- tryCatch(names(raw_df()), error = function(e) character(0))

  default_peak <- if (!is.null(input$row_id_col) && input$row_id_col %in% raw_cols) {
    input$row_id_col
  } else {
    guess_col_ci(
      raw_cols,
      c("row ID", "row id", "id", "feature_id"),
      default = "__none__"
    )
  }

  tagList(
    selectInput(
      "gnps_cluster1_col",
      "GNPS ClusterID1 column:",
      choices = g_cols,
      selected = guess_col_ci(g_cols, c("ClusterID1", "CLUSTERID1"), g_cols[1])
    ),

    selectInput(
      "gnps_cluster2_col",
      "GNPS ClusterID2 column:",
      choices = g_cols,
      selected = guess_col_ci(g_cols, c("ClusterID2", "CLUSTERID2"), g_cols[min(2, length(g_cols))])
    ),

    selectInput(
      "gnps_component_col",
      "GNPS ComponentIndex column:",
      choices = g_cols,
      selected = guess_col_ci(g_cols, c("ComponentIndex", "component"), g_cols[1])
    ),

    selectInput(
      "gnps_peak_id_col",
      "Optional peak-table pruning by Peak ID:",
      choices = c("Do not use peak-table ID matching" = "__none__", raw_cols),
      selected = "__none__"
    )
  )
})

sirius_stats_data <- reactive({
  req(sirius_df(), input$stats_sirius_id_col, input$stats_sirius_col)

  s <- as.data.frame(
    sirius_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  validate(
    need(input$stats_sirius_id_col %in% names(s), "Selected SIRIUS ID column was not found."),
    need(input$stats_sirius_col %in% names(s), "Selected SIRIUS statistics column was not found.")
  )

  out <- tibble::tibble(
    SIRIUS_ID = trimws(as.character(s[[input$stats_sirius_id_col]])),
    Annotation = clean_stats_value(s[[input$stats_sirius_col]])
  ) %>%
    dplyr::filter(nzchar(SIRIUS_ID), !is.na(Annotation))

  keep_ids <- selected_peak_ids_for_gnps()

  if (!is.null(keep_ids) && length(keep_ids)) {
    out <- out %>%
      dplyr::filter(SIRIUS_ID %in% keep_ids)
  }

  out
})

output$stats_class_picker <- renderUI({
  req(sirius_stats_data())

  vals <- sirius_stats_data() %>%
    dplyr::count(Annotation, sort = TRUE, name = "Frequency") %>%
    dplyr::pull(Annotation)

  if (!length(vals)) {
    return(
      div(
        class = "small-note",
        "No non-empty values detected for the selected SIRIUS column."
      )
    )
  }

  pickerInput(
    "stats_selected_class",
    "Specific class/value for GNPS ComponentIndex statistics:",
    choices = vals,
    selected = vals[1],
    multiple = FALSE,
    options = list(
      `live-search` = TRUE,
      `style` = "btn-success"
    )
  )
})

sirius_frequency_table <- reactive({
  req(sirius_stats_data())

  n_total <- nrow(sirius_stats_data())

  sirius_stats_data() %>%
    dplyr::group_by(Annotation) %>%
    dplyr::summarise(
      Frequency = dplyr::n(),
      Unique_IDs = dplyr::n_distinct(SIRIUS_ID),
      Percent = round(100 * Frequency / n_total, 2),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Frequency), Annotation)
})

output$sirius_frequency_table <- DT::renderDT({
  datatable(
    sirius_frequency_table(),
    rownames = FALSE,
    class = "compact stripe hover nowrap",
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      order = list(list(1, "desc"))
    )
  )
})

selected_peak_ids_for_gnps <- reactive({
  if (is.null(input$gnps_peak_id_col) || identical(input$gnps_peak_id_col, "__none__")) {
    return(NULL)
  }

  req(raw_df())

  raw <- as.data.frame(
    raw_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  validate(
    need(input$gnps_peak_id_col %in% names(raw), "Selected peak-table ID column was not found.")
  )

  ids <- trimws(as.character(raw[[input$gnps_peak_id_col]]))
  ids[nzchar(ids)]
})

gnps_component_map <- reactive({
  req(
    gnps_pairs_df(),
    input$gnps_cluster1_col,
    input$gnps_cluster2_col,
    input$gnps_component_col
  )

  g <- as.data.frame(
    gnps_pairs_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  validate(
    need(input$gnps_cluster1_col %in% names(g), "GNPS ClusterID1 column was not found."),
    need(input$gnps_cluster2_col %in% names(g), "GNPS ClusterID2 column was not found."),
    need(input$gnps_component_col %in% names(g), "GNPS ComponentIndex column was not found.")
  )

  out <- g %>%
    dplyr::transmute(
      ComponentIndex = as.character(.data[[input$gnps_component_col]]),
      ClusterID1 = as.character(.data[[input$gnps_cluster1_col]]),
      ClusterID2 = as.character(.data[[input$gnps_cluster2_col]])
    ) %>%
    tidyr::pivot_longer(
      cols = c("ClusterID1", "ClusterID2"),
      names_to = "Cluster_side",
      values_to = "ClusterID"
    ) %>%
    dplyr::mutate(
      ClusterID = trimws(as.character(ClusterID)),
      ComponentIndex = trimws(as.character(ComponentIndex))
    ) %>%
    dplyr::filter(nzchar(ClusterID), nzchar(ComponentIndex)) %>%
    dplyr::distinct(ComponentIndex, ClusterID)

  keep_ids <- selected_peak_ids_for_gnps()

  if (!is.null(keep_ids) && length(keep_ids)) {
    out <- out %>%
      dplyr::filter(ClusterID %in% keep_ids)
  }

  out
})

selected_class_component_stats <- reactive({
  req(sirius_stats_data(), input$stats_selected_class)

  ids <- sirius_stats_data() %>%
    dplyr::filter(Annotation == input$stats_selected_class) %>%
    dplyr::distinct(SIRIUS_ID)

  if (!nrow(ids)) {
    return(
      tibble::tibble(
        ComponentIndex = character(),
        Points = integer(),
        ClusterIDs = character()
      )
    )
  }

  if (is.null(input$file_gnps_pairs)) {
    return(
      tibble::tibble(
        ComponentIndex = "GNPS pairs not uploaded",
        Points = dplyr::n_distinct(ids$SIRIUS_ID),
        ClusterIDs = paste(sort(unique(ids$SIRIUS_ID)), collapse = ", ")
      )
    )
  }

  comp <- gnps_component_map()

  ids %>%
    dplyr::left_join(comp, by = c("SIRIUS_ID" = "ClusterID")) %>%
    dplyr::mutate(
      ComponentIndex = dplyr::if_else(
        is.na(ComponentIndex) | !nzchar(ComponentIndex),
        "No ComponentIndex match",
        ComponentIndex
      )
    ) %>%
    dplyr::group_by(ComponentIndex) %>%
    dplyr::summarise(
      Points = dplyr::n_distinct(SIRIUS_ID),
      ClusterIDs = paste(sort(unique(SIRIUS_ID)), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Points), ComponentIndex)
})

output$selected_class_summary <- renderUI({
  req(sirius_stats_data(), input$stats_selected_class)

  ids_all <- sirius_stats_data() %>%
    dplyr::filter(Annotation == input$stats_selected_class) %>%
    dplyr::distinct(SIRIUS_ID)

  comp_stats <- selected_class_component_stats()

  mapped <- if (!is.null(input$file_gnps_pairs)) {
    sum(comp_stats$Points[comp_stats$ComponentIndex != "No ComponentIndex match"], na.rm = TRUE)
  } else {
    NA_integer_
  }

  div(
    style = "
      background:#ffffffcc;
      padding:10px;
      border-radius:10px;
      border:1px solid #ddd;
      margin-bottom:10px;
    ",
    h4("Selected class summary"),
    tags$ul(
      tags$li(strong("Column: "), input$stats_sirius_col),
      tags$li(strong("Selected value: "), input$stats_selected_class),
      tags$li(strong("Unique SIRIUS IDs / points: "), dplyr::n_distinct(ids_all$SIRIUS_ID)),
      tags$li(
        strong("SIRIUS entries: "),
        nrow(
          sirius_stats_data() %>%
            dplyr::filter(Annotation == input$stats_selected_class)
        )
      ),
      if (!is.na(mapped)) {
        tags$li(strong("Points with ComponentIndex match: "), mapped)
      },
      if (!is.null(input$file_gnps_pairs)) {
        tags$li(strong("Number of ComponentIndex groups shown: "), nrow(comp_stats))
      }
    )
  )
})

output$selected_class_component_table <- DT::renderDT({
  datatable(
    selected_class_component_stats(),
    rownames = FALSE,
    class = "compact stripe hover nowrap",
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      order = list(list(1, "desc"))
    )
  )
})

output$sirius_gnps_main <- renderUI({
  if (!isTRUE(input$use_sirius) || is.null(input$file_sirius)) {
    return(
      div(
        class = "highlight",
        "No SIRIUS annotation table uploaded yet. Go to Load & Process -> Join with Annotation -> upload SIRIUS summary."
      )
    )
  }

  tagList(
    h3("Frequency by selected SIRIUS column"),
    withSpinner(DTOutput("sirius_frequency_table"), type = 8, color = "#66CDAA"),
    tags$hr(),
    uiOutput("selected_class_summary"),
    h3("Selected class distribution by GNPS ComponentIndex"),
    withSpinner(DTOutput("selected_class_component_table"), type = 8, color = "#66CDAA")
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
      scale_data = isTRUE(input$standard_scaling),
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
      updateMaterialSwitch(session, "use_npc_filter", value = FALSE)
      updateMaterialSwitch(session, "use_classyfire_filter", value = FALSE)

      # pickers / radios
      updatePickerInput(session, "sel_feat", selected = character(0))
      updatePickerInput(session, "npc_filter_values", selected = character(0))
      updatePickerInput(session, "classyfire_filter_values", selected = character(0))
      updateRadioButtons(session, "color_by", selected = "Groups")
      updateRadioButtons(session, "present_as", selected = "Boxplot")
      updateRadioButtons(session, "fc_dir", selected = "both")

      # sliders
      dd <- rv$volcano
      mzr <- finite_range(dd$mz)
      rtr <- finite_range(dd$RT)
      mr  <- finite_range(log10(dd$Mean + 1.1))

      if (!is.null(mzr)) {
        updateSliderInput(
          session, "mz_range",
          value = c(
            floor(mzr[1] * 10000) / 10000,
            ceiling(mzr[2] * 10000) / 10000
          )
        )
      }

      if (!is.null(rtr)) {
        updateSliderInput(
          session, "rt_range",
          value = c(
            floor(rtr[1] * 1000) / 1000,
            ceiling(rtr[2] * 1000) / 1000
          )
        )
      }

      if (!is.null(mr)) {
        updateSliderInput(
          session, "intensity_range",
          value = c(
            floor(mr[1] * 1000) / 1000,
            ceiling(mr[2] * 1000) / 1000
          )
        )
      }

      updateNumericInput(session, "sig_p_cutoff", value = 0.05)
      updateSliderInput(session, "fc_thr", value = 1)
      updateSelectInput(session, "volcano_palette", selected = "Set1")
      updateSelectInput(session, "box_palette", selected = "Dark2")
    })

  annotation_filter_choices <- reactive({
  req(procReady())

  dd <- rv$volcano

  clean_choices <- function(x) {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x) & x != "Not provided"]
    sort(unique(x))
  }

  list(
    npc = clean_choices(dd$`NPC#class`),
    classyfire = clean_choices(dd$`ClassyFire#class`)
  )
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
      materialSwitch("sig_only", "Significant only", value = FALSE, status = "success"),

tags$hr(),
h4(class = "highlight", "Annotation filters"),

materialSwitch(
  "use_npc_filter",
  "Filter by NPC class",
  value = FALSE,
  status = "success"
),

conditionalPanel(
  condition = "input.use_npc_filter == true",
  if (length(annotation_filter_choices()$npc) > 0) {
    pickerInput(
      inputId = "npc_filter_values",
      label = "NPC class(es):",
      choices = annotation_filter_choices()$npc,
      selected = character(0),
      multiple = TRUE,
      options = list(
  `actions-box` = TRUE,
  `live-search` = TRUE,
  `none-selected-text` = "Select NPC class(es)",
  `style` = "btn-success",
  `selected-text-format` = "count > 1",
  `count-selected-text` = "{0} NPC class(es) selected"
)
    )
  } else {
    div(
      class = "small-note",
      "No NPC classes detected. Upload SIRIUS annotation and rerun preprocessing."
    )
  }
),

materialSwitch(
  "use_classyfire_filter",
  "Filter by ClassyFire class",
  value = FALSE,
  status = "success"
),

conditionalPanel(
  condition = "input.use_classyfire_filter == true",
  if (length(annotation_filter_choices()$classyfire) > 0) {
    pickerInput(
      inputId = "classyfire_filter_values",
      label = "ClassyFire class(es):",
      choices = annotation_filter_choices()$classyfire,
      selected = character(0),
      multiple = TRUE,
      options = list(
  `actions-box` = TRUE,
  `live-search` = TRUE,
  `none-selected-text` = "Select ClassyFire class(es)",
  `style` = "btn-success",
  `selected-text-format` = "count > 1",
  `count-selected-text` = "{0} ClassyFire class(es) selected"
)
    )
  } else {
    div(
      class = "small-note",
      "No ClassyFire classes detected. Upload SIRIUS annotation and rerun preprocessing."
    )
  }
),

tags$hr(),

pickerInput(
  inputId = "sel_feat",
  label   = "Select/deselect features (optional):",
  choices = sort(unique(rv$volcano$Feature)),
  options = list(
  `actions-box` = TRUE,
  `live-search` = TRUE,
  `style` = "btn-success",
  `selected-text-format` = "count > 2",
  `count-selected-text` = "{0} feature(s) selected"
),
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

selectInput(
  "volcano_palette",
  "Volcano palette:",
  choices = palette_choices,
  selected = "Set1"
),

selectInput(
  "box_palette",
  "Boxplot/scatter palette:",
  choices = palette_choices,
  selected = "Dark2"
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
  mr  <- finite_range(log10(dd$Mean + 1.1))
  fcmax <- max(abs(dd$FC), na.rm = TRUE)

  validate(need(!is.null(mzr) && !is.null(rtr) && !is.null(mr),
                "No finite mz/RT/Mean values available for sliders."))

  mz_min <- floor(mzr[1] * 10000) / 10000
  mz_max <- ceiling(mzr[2] * 10000) / 10000

  rt_min <- floor(rtr[1] * 1000) / 1000
  rt_max <- ceiling(rtr[2] * 1000) / 1000

  int_min <- floor(mr[1] * 1000) / 1000
  int_max <- ceiling(mr[2] * 1000) / 1000

  tagList(
    sliderInput("mz_range", "m/z:",
                min = mz_min,
                max = mz_max,
                value = c(mz_min, mz_max),
                step = 0.0001),

    sliderInput("rt_range", "RT:",
                min = rt_min,
                max = rt_max,
                value = c(rt_min, rt_max),
                step = 0.001),

    sliderInput("intensity_range", "Mean log10(Intensity):",
                min = int_min,
                max = int_max,
                value = c(int_min, int_max),
                step = 0.001),

          tags$hr(),

h4(class = "highlight", "Significance thresholds"),

numericInput(
  "sig_p_cutoff",
  "Adj.p-value cutoff:",
  value = 0.05,
  min = 0,
  max = 1,
  step = 0.001
),

sliderInput(
  "fc_thr",
  "FC threshold (|log2FC| ≥):",
  min = 0,
  max = max(1, round(fcmax, 1)),
  value = 1,
  step = 0.1
),

materialSwitch(
  "use_fdr_filter",
  "Filter by Adj.p-value threshold",
  value = FALSE,
  status = "success"
),

materialSwitch(
  "use_fc_filter",
  "Filter by Fold-Change threshold",
  value = FALSE,
  status = "success"
),

conditionalPanel(
  condition = "input.use_fc_filter == true",
  radioButtons(
    "fc_dir",
    "Direction:",
    choices = c("Both sides" = "both", "Up only" = "up", "Down only" = "down"),
    selected = "both",
    inline = TRUE
  )
)
    )
  })

  # ---- Filtered volcano data
  filtered_volcano <- reactive({
    req(procReady(), input$mz_range, input$rt_range, input$intensity_range)

    dd <- rv$volcano
    pcut <- suppressWarnings(as.numeric(input$sig_p_cutoff %||% 0.05))
    if (!is.finite(pcut) || pcut <= 0 || pcut > 1) pcut <- 0.05

    fcut <- suppressWarnings(as.numeric(input$fc_thr %||% 1))
    if (!is.finite(fcut) || fcut < 0) fcut <- 1

    if (!is.null(input$sel_feat) && length(input$sel_feat) > 0) {
  dd <- dd %>% dplyr::filter(Feature %in% input$sel_feat)
}
    if (isTRUE(input$sig_only)) {
  dd <- dd %>%
    dplyr::filter(
      is.finite(`Adj.p-value`),
      `Adj.p-value` <= pcut,
      is.finite(FC),
      abs(FC) >= fcut
    )
}

    if (isTRUE(input$use_npc_filter)) {
  validate(
    need(
      !is.null(input$npc_filter_values) && length(input$npc_filter_values) > 0,
      "NPC filter is enabled. Select at least one NPC class."
    )
  )

  dd <- dd %>%
    dplyr::filter(`NPC#class` %in% input$npc_filter_values)
}

if (isTRUE(input$use_classyfire_filter)) {
  validate(
    need(
      !is.null(input$classyfire_filter_values) && length(input$classyfire_filter_values) > 0,
      "ClassyFire filter is enabled. Select at least one ClassyFire class."
    )
  )

  dd <- dd %>%
    dplyr::filter(`ClassyFire#class` %in% input$classyfire_filter_values)
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
  dd <- dd %>%
    dplyr::filter(
      is.finite(`Adj.p-value`),
      `Adj.p-value` <= pcut
    )
}

    if (isTRUE(input$use_fc_filter)) {
      thr <- fcut
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

    fc_line <- suppressWarnings(as.numeric(input$fc_thr %||% 1))
    if (!is.finite(fc_line) || fc_line < 0) fc_line <- 1

    p_thr <- suppressWarnings(as.numeric(input$sig_p_cutoff %||% 0.05))
    if (!is.finite(p_thr) || p_thr <= 0 || p_thr > 1) p_thr <- 0.05

    ythr <- -log10(pmax(p_thr, .Machine$double.xmin))

    shapes <- list(
  list(
    type = "line",
    x0 = -fc_line,
    x1 = -fc_line,
    xref = "x",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    line = list(dash = "dot")
  ),
  list(
    type = "line",
    x0 = fc_line,
    x1 = fc_line,
    xref = "x",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    line = list(dash = "dot")
  ),
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = ythr,
    y1 = ythr,
    yref = "y",
    line = list(dash = "dot")
  )
)

    if (input$color_by == "Groups") {
      plot_ly(
        data = dd,
        x = ~FC, y = ~`Adj.p-value.log`,
        colors = make_palette(input$volcano_palette %||% "Set1", length(unique(dd$Groups))),
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
        colors = make_palette(input$volcano_palette %||% "viridis", 12),
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

    box_cols <- make_palette(input$box_palette %||% "Dark2", length(unique(xx)))

    if (input$present_as == "Boxplot") {
      plot_ly(
        x = xx, y = yy,
        colors = box_cols,
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
        colors = box_cols,
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
  filename = function() {
    req(rv$volcano)

    n_comp <- length(unique(rv$volcano$Groups))

    if (n_comp > 1) {
      paste0(dataset_name(), "_volcano_table_wide.csv")
    } else {
      paste0(dataset_name(), "_volcano_table.csv")
    }
  },

  content = function(file) {
    req(rv$volcano)

    out <- volcano_to_wide_if_needed(rv$volcano)

    data.table::fwrite(out, file, na = "")
  }
)

output$dl_matrix <- downloadHandler(
  filename = function() {
    paste0(dataset_name(), "_MetaboAnalyst_table.csv")
  },
  content = function(file) {
    req(rv$df_used, rv$mat)

    out <- as.data.frame(rv$df_used, check.names = FALSE, stringsAsFactors = FALSE)
    out <- cbind(
      Sample = rownames(rv$mat),
      out
    )

    data.table::fwrite(out, file, na = "")
  }
)

output$dl_autoplotter_zip <- downloadHandler(
  filename = function() {
    nm <- input$file_data$name %||% "dataset.csv"
    paste0(tools::file_path_sans_ext(basename(nm)), "_AutoPlotter.zip")
  },

  content = function(file) {
    req(rv$df_used, rv$mat, rv$fmap)

    zip_dir <- tempfile("autoplotter_")
    dir.create(zip_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(zip_dir, recursive = TRUE, force = TRUE), add = TRUE)

    data_file <- file.path(zip_dir, "data_table.csv")
    meta_file <- file.path(zip_dir, "metadata.csv")
    name_file <- file.path(zip_dir, "name_map.csv")

    autoplotter_data <- make_autoplotter_data(
      df_used = rv$df_used,
      sample_names = rownames(rv$mat)
    )

    autoplotter_metadata <- make_autoplotter_metadata(
      df_used = rv$df_used,
      sample_names = rownames(rv$mat)
    )

    autoplotter_name_map <- make_autoplotter_name_map(
      fmap = rv$fmap,
      volcano = rv$volcano
    )

    data.table::fwrite(autoplotter_data, data_file, na = "")
    data.table::fwrite(autoplotter_metadata, meta_file, na = "")
    data.table::fwrite(autoplotter_name_map, name_file, na = "")

    tmp_zip <- tempfile(fileext = ".zip")
    on.exit(unlink(tmp_zip, force = TRUE), add = TRUE)

    zip::zipr(
      zipfile = tmp_zip,
      files = c(data_file, meta_file, name_file),
      root = zip_dir
    )

    ok <- file.copy(tmp_zip, file, overwrite = TRUE)
    if (!ok) {
      stop("Failed to copy AutoPlotter ZIP archive to download file.")
    }
  },

  contentType = "application/zip"
)

}
