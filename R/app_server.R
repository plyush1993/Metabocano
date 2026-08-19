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

  observeEvent(
  input$selected_feature_info_copied,
  {
    showNotification(
      "m/z and RT copied as: mz,rt",
      type = "message",
      duration = 2
    )
  },
  ignoreInit = TRUE
)

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

  output$manual_sample_cols_ui <- renderUI({

  req(raw_df())

  selectizeInput(
    "sample_cols_manual",
    "Pick sample columns:",
    choices = names(raw_df()),
    selected = NULL,
    multiple = TRUE,
    options = list(
      placeholder = "Select first and last sample columns"
    )
  )
})

  observeEvent(input$sample_cols_manual, {

  req(raw_df())

  sel <- input$sample_cols_manual

  if (is.null(sel) || length(sel) < 2)
    return()

  cols <- names(raw_df())

  pos <- match(sel, cols)
  pos <- pos[!is.na(pos)]

  if (length(pos) < 2)
    return()

  # Select everything between leftmost and rightmost selection
  range_cols <- cols[min(pos):max(pos)]

  if (!setequal(sel, range_cols)) {

    updateSelectizeInput(
      session,
      "sample_cols_manual",
      selected = range_cols
    )
  }

}, ignoreInit = TRUE)

  output$raw_header <- renderUI({
    req(raw_df())
    h3(sprintf("Raw dataset: %d rows × %d columns", nrow(raw_df()), ncol(raw_df())))
  })

  output$raw_preview <- renderDT({
    req(raw_df())
    datatable(head(raw_df(), 20), options = list(scrollX = TRUE, pageLength = 8))
  })

  sample_cols_selected <- reactive({

  req(raw_df(), input$row_id_col, input$mz_col, input$rt_col)

  df <- as.data.frame(
    raw_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cols <- names(df)

  mode <- input$sample_mode %||% "kws"

  # Columns that must never be treated as samples in Auto/Keyword mode
  meta <- unique(c(
    input$row_id_col,
    input$mz_col,
    input$rt_col,
    "<Auto-generate>"
  ))

  meta <- meta[
    !is.na(meta) &
      nzchar(meta)
  ]

  # MANUAL
  if (mode == "manual") {

    validate(
      need(
        !is.null(input$sample_cols_manual) &&
          length(input$sample_cols_manual) > 0,
        "Pick sample columns."
      )
    )

    sc <- intersect(
      input$sample_cols_manual,
      cols
    )

    validate(
      need(
        length(sc) > 0,
        "Manual sample columns not found in table."
      )
    )

    return(sc)
  }

  # KEYWORDS
  if (mode == "kws") {

    kws <- input$sample_keywords %||%
      character(0)

    idx <- multi_sample_idx(
      cols,
      kws
    )

    validate(
      need(
        length(idx) > 0,
        paste0(
          "No sample columns matched the keywords: ",
          paste(kws, collapse = ", ")
        )
      )
    )

    sc <- cols[idx]

    # Do not allow m/z, RT, ID etc. to become samples
    sc <- setdiff(
      sc,
      meta
    )

    validate(
      need(
        length(sc) > 0,
        "Keyword hits were only metadata columns. Use Manual or Auto."
      )
    )

    return(sc)
  }

  # AUTO
  cand <- setdiff(
    cols,
    meta
  )

  # Same idea as MetaboCensoR:
  # exclude typical 'row ...' helper columns
  cand <- cand[
    !grepl(
      "^row\\b",
      cand,
      ignore.case = TRUE
    )
  ]

  prop_num <- vapply(
    df[cand],
    function(x) {

      x2 <- suppressWarnings(
        as.numeric(
          as.character(x)
        )
      )

      mean(
        is.finite(x2),
        na.rm = TRUE
      )
    },
    numeric(1)
  )

  # At least 70% numeric-like
  sc <- cand[
    prop_num >= 0.7
  ]

  validate(
    need(
      length(sc) > 0,
      "Auto-detect found no numeric sample columns. Switch to Manual or Keywords."
    )
  )

  sc
})

  output$peak_extra_cols_ui <- renderUI({

  req(raw_df())

  df <- raw_df()
  cols <- names(df)

  # Identify sample-intensity columns so they are not offered
  # as additional feature metadata columns.
  sample_cols <- sample_cols_selected()

  # Exclude columns already used as essential feature metadata
  excluded_cols <- unique(
    c(
      input$row_id_col %||% character(0),
      input$mz_col %||% character(0),
      input$rt_col %||% character(0),
      sample_cols,
      "<Auto-generate>"
    )
  )

  choices <- setdiff(
    cols,
    excluded_cols
  )

  if (!length(choices)) {

    return(
      div(
        class = "small-note",
        "No additional non-sample peak-table columns were detected."
      )
    )
  }

  current_selection <- isolate(
    input$peak_extra_cols
  ) %||% character(0)

  current_selection <- intersect(
    current_selection,
    choices
  )

  pickerInput(
    inputId = "peak_extra_cols",
    label = "Peak-table columns to include:",
    choices = choices,
    selected = current_selection,
    multiple = TRUE,

    options = list(
      `actions-box` = TRUE,
      `live-search` = TRUE,
      `none-selected-text` =
        "Select one or more peak-table columns",
      `selected-text-format` = "count > 2",
      `count-selected-text` =
        "{0} peak-table column(s) selected",
      `style` = "btn-success"
    )
  )
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

  sc <- sample_cols_selected()

  parse_feature_table_to_matrix(
    df,
    row_id_col = rid_col,
    mz_col = input$mz_col,
    rt_col = input$rt_col,
    sample_cols = sc
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

  labels_from_sample_names_or_raw(
    sample_names(),
    token_sep = input$token_sep %||% "_",
    token_index = input$token_index %||% 2,
    clean_names = FALSE
  )
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

  comparison_pairs <- reactive({
  req(labels_vec())

  levs <- sort(
    unique(
      trimws(as.character(labels_vec()))
    )
  )

  levs <- levs[nzchar(levs)]

  validate(
    need(
      length(levs) >= 2,
      "Need at least 2 groups to define comparisons."
    )
  )

  # All directional comparisons:
  # A / B and B / A are separate options
  tidyr::expand_grid(
    Group_num = levs,
    Group_den = levs
  ) %>%
    dplyr::filter(Group_num != Group_den) %>%
    dplyr::mutate(
      Comparison_ID = sprintf(
        "comparison_%04d",
        dplyr::row_number()
      ),
      Comparison = paste0(
        Group_num,
        " / ",
        Group_den
      )
    ) %>%
    dplyr::select(
      Comparison_ID,
      Comparison,
      Group_num,
      Group_den
    )
})


output$comparison_picker <- renderUI({
  req(labels_vec())

  levs <- sort(
    unique(
      trimws(as.character(labels_vec()))
    )
  )

  pairs <- comparison_pairs()

  # Preserve the currently selected reference group
  old_ref <- isolate(input$ref_group)

  if (is.null(old_ref) || !old_ref %in% levs) {
    old_ref <- levs[1]
  }

  # Preserve valid manual choices when UI is rebuilt
  old_manual <- isolate(
    input$manual_comparisons
  ) %||% character(0)

  old_manual <- intersect(
    old_manual,
    pairs$Comparison_ID
  )

  tagList(

    conditionalPanel(
      condition = "input.comparison_mode == 'reference'",

      selectInput(
        "ref_group",
        "Reference Group:",
        choices = levs,
        selected = old_ref
      )
    ),

    conditionalPanel(
      condition = "input.comparison_mode == 'manual'",

      pickerInput(
        inputId = "manual_comparisons",
        label = "Select comparisons:",
        choices = stats::setNames(
          pairs$Comparison_ID,
          pairs$Comparison
        ),
        selected = old_manual,
        multiple = TRUE,
        options = list(
          `actions-box` = TRUE,
          `live-search` = TRUE,
          `none-selected-text` =
            "Select one or more comparisons",
          `selected-text-format` = "count > 2",
          `count-selected-text` =
            "{0} comparison(s) selected",
          `style` = "btn-success"
        )
      )
    )
  )
})

selected_manual_comparisons <- reactive({
  req(
    identical(
      input$comparison_mode %||% "reference",
      "manual"
    )
  )

  selected_ids <- input$manual_comparisons %||%
    character(0)

  validate(
    need(
      length(selected_ids) > 0,
      "Select at least one manual comparison."
    )
  )

  pairs <- comparison_pairs()

  selected_rows <- match(
    selected_ids,
    pairs$Comparison_ID
  )

  validate(
    need(
      !any(is.na(selected_rows)),
      paste0(
        "One or more selected comparisons ",
        "are no longer available."
      )
    )
  )

  pairs[
    selected_rows,
    c("Group_num", "Group_den"),
    drop = FALSE
  ]
})

  # ---- SIRIUS pickers ----
  sirius_df <- reactive({

  req(input$use_sirius)
  req(input$file_sirius)

  ext <- tolower(
    tools::file_ext(
      input$file_sirius$name
    )
  )

  validate(
    need(
      ext %in% c("csv", "tsv", "txt"),
      "SIRIUS file must be .csv, .tsv, or .txt."
    )
  )

  delim <- if (identical(ext, "csv")) {
    ","
  } else {
    "\t"
  }

  as.data.frame(
    vroom::vroom(
      input$file_sirius$datapath,
      delim = delim,
      col_names = TRUE,
      show_col_types = FALSE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})

  output$sirius_pickers <- renderUI({

  req(sirius_df())

  cols <- names(
    sirius_df()
  )

  tagList(

    selectInput(
      "sirius_idcol",
      "SIRIUS Feature ID column:",
      choices = cols,
      selected = if (
        "mappingFeatureId" %in% cols
      ) {
        "mappingFeatureId"
      } else {
        cols[1]
      }
    ),

    selectInput(
      "sirius_npcol",
      "NPC column:",
      choices = cols,
      selected = if (
        "NPC#class" %in% cols
      ) {
        "NPC#class"
      } else {
        cols[1]
      }
    ),

    selectInput(
      "sirius_cfcol",
      "ClassyFire column:",
      choices = cols,
      selected = if (
        "ClassyFire#class" %in% cols
      ) {
        "ClassyFire#class"
      } else {
        cols[1]
      }
    ),

    materialSwitch(
      inputId = "use_sirius_extra_cols",
      label = "Add additional SIRIUS columns",
      value = FALSE,
      status = "success",
      width = "auto"
    ),

    conditionalPanel(
      condition = "input.use_sirius_extra_cols == true",

      uiOutput(
        "sirius_extra_cols_ui"
      )
    )
  )
})

  output$sirius_extra_cols_ui <- renderUI({

  req(
    sirius_df(),
    input$sirius_idcol,
    input$sirius_npcol,
    input$sirius_cfcol
  )

  cols <- names(
    sirius_df()
  )

  choices <- setdiff(
    cols,
    c(
      input$sirius_idcol,
      input$sirius_npcol,
      input$sirius_cfcol
    )
  )

  if (!length(choices)) {

    return(
      div(
        class = "small-note",
        "No additional SIRIUS columns are available."
      )
    )
  }

  current_selection <- isolate(
    input$sirius_extra_cols
  ) %||% character(0)

  current_selection <- intersect(
    current_selection,
    choices
  )

  pickerInput(
    inputId = "sirius_extra_cols",
    label = "Additional SIRIUS columns:",
    choices = choices,
    selected = current_selection,
    multiple = TRUE,

    options = list(
      `actions-box` = TRUE,
      `live-search` = TRUE,
      `none-selected-text` =
        "Select one or more SIRIUS columns",
      `selected-text-format` = "count > 2",
      `count-selected-text` =
        "{0} SIRIUS column(s) selected",
      `style` = "btn-success"
    )
  )
})

gnps_annotation_df <- reactive({

  req(input$use_gnps_annotation)
  req(input$file_gnps_annotation)

  ext <- tolower(
    tools::file_ext(
      input$file_gnps_annotation$name
    )
  )

  validate(
    need(
      ext %in% c("tsv", "txt", "csv"),
      "GNPS annotation file must be .tsv, .txt, or .csv."
    )
  )

  delim <- if (identical(ext, "csv")) {
    ","
  } else {
    "\t"
  }

  as.data.frame(
    vroom::vroom(
      input$file_gnps_annotation$datapath,
      delim = delim,
      col_names = TRUE,
      show_col_types = FALSE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})

output$gnps_annotation_pickers <- renderUI({

  req(gnps_annotation_df())

  cols <- names(
    gnps_annotation_df()
  )

  validate(
    need(
      length(cols) > 0,
      "No columns were detected in the GNPS file."
    )
  )

  # Default GNPS ID column
  default_id <- guess_col(
    cols,
    c(
      "#Scan#",
      "Scan",
      "scan",
      "ClusterIndex",
      "Cluster ID",
      "row ID",
      "id"
    )
  ) %||% cols[1]

  # Default primary annotation column
  default_annotation <- guess_col(
    cols,
    c(
      "Compound_Name",
      "Compound_name",
      "Compound name",
      "CompoundName",
      "Library compound name",
      "Annotation",
      "Name"
    )
  ) %||% cols[1]

  tagList(

    selectInput(
      "gnps_annotation_idcol",
      "GNPS ID column:",
      choices = cols,
      selected = default_id
    ),

    selectInput(
      "gnps_annotation_col",
      "Primary GNPS annotation column:",
      choices = cols,
      selected = default_annotation
    ),

    materialSwitch(
      inputId = "use_gnps_extra_cols",
      label = "Add additional GNPS columns",
      value = FALSE,
      status = "success",
      width = "auto"
    ),

    conditionalPanel(
      condition = "input.use_gnps_extra_cols == true",

      uiOutput(
        "gnps_extra_cols_ui"
      )
    )
  )
})

output$gnps_extra_cols_ui <- renderUI({

  req(
    gnps_annotation_df(),
    input$gnps_annotation_idcol,
    input$gnps_annotation_col
  )

  cols <- names(
    gnps_annotation_df()
  )

  # Do not offer the join ID or primary annotation again
  choices <- setdiff(
    cols,
    c(
      input$gnps_annotation_idcol,
      input$gnps_annotation_col
    )
  )

  if (!length(choices)) {

    return(
      div(
        class = "small-note",
        "No additional GNPS columns are available."
      )
    )
  }

  current_selection <- isolate(
    input$gnps_extra_cols
  ) %||% character(0)

  current_selection <- intersect(
    current_selection,
    choices
  )

  pickerInput(
    inputId = "gnps_extra_cols",
    label = "Additional GNPS columns:",
    choices = choices,
    selected = current_selection,
    multiple = TRUE,

    options = list(
      `actions-box` = TRUE,
      `live-search` = TRUE,
      `none-selected-text` =
        "Select one or more GNPS columns",
      `selected-text-format` = "count > 2",
      `count-selected-text` =
        "{0} GNPS column(s) selected",
      `style` = "btn-success"
    )
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

materialSwitch(
  inputId = "prune_sirius_to_peak",
  label = "Restrict SIRIUS statistics to uploaded peak table",
  value = TRUE,
  status = "success",
  width = "auto"
),

uiOutput(
  "sirius_pruning_notice"
),

uiOutput(
  "sirius_peak_id_picker"
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
    uiOutput("network_component_picker"),

    div(
      class = "small-note",
      "GNPS pairs are converted from ClusterID1/ClusterID2 into one ClusterID column. By default, no peak-table pruning is applied. If a peak-table ID column is selected, SIRIUS IDs and GNPS ClusterIDs are restricted to IDs present in the uploaded peak table."
    )
  )
})

output$sirius_pruning_notice <- renderUI({

  if (!isTRUE(input$prune_sirius_to_peak)) {
    return(NULL)
  }

  if (!is.null(input$file_data)) {
    return(NULL)
  }

})

output$sirius_peak_id_picker <- renderUI({

    if (!isTRUE(input$prune_sirius_to_peak)) {
    return(NULL)
  }

  if (is.null(input$file_data)) {
    return(NULL)
  }

  req(raw_df())

  raw_cols <- names(
    raw_df()
  )

  validate(
    need(
      length(raw_cols) > 0,
      "No peak-table columns were detected."
    )
  )

  selected_row_id <- input$row_id_col %||%
    ""

  default_peak_id <- if (
    identical(
      selected_row_id,
      "<Auto-generate>"
    )
  ) {
    "__auto__"
  } else if (
    selected_row_id %in% raw_cols
  ) {
    selected_row_id
  } else {
    guess_col_ci(
      raw_cols,
      c(
        "row ID",
        "row id",
        "id",
        "feature_id",
        "feature id"
      ),
      default = raw_cols[1]
    )
  }

pickerInput(
  inputId = "sirius_peak_id_col",
  label = "Peak-table column used for SIRIUS pruning:",

  choices = c(
    "Auto-generated feature ID (row number)" = "__auto__",
    raw_cols
  ),

  selected = default_peak_id,
  multiple = FALSE,

  options = list(
    `live-search` = TRUE,
    `size` = 10,
    `style` = "btn-success",
    `none-selected-text` = "Choose a peak-table ID column"
  )
)
})

output$gnps_pickers_ui <- renderUI({

  req(
    gnps_pairs_df()
  )

  g_cols <- names(
    gnps_pairs_df()
  )

  tagList(

    selectInput(
      "gnps_cluster1_col",
      "GNPS ClusterID1 column:",
      choices = g_cols,
      selected = guess_col_ci(
        g_cols,
        c(
          "ClusterID1",
          "CLUSTERID1"
        ),
        g_cols[1]
      )
    ),

    selectInput(
      "gnps_cluster2_col",
      "GNPS ClusterID2 column:",
      choices = g_cols,
      selected = guess_col_ci(
        g_cols,
        c(
          "ClusterID2",
          "CLUSTERID2"
        ),
        g_cols[
          min(
            2,
            length(g_cols)
          )
        ]
      )
    ),

    selectInput(
      "gnps_component_col",
      "GNPS ComponentIndex column:",
      choices = g_cols,
      selected = guess_col_ci(
        g_cols,
        c(
          "ComponentIndex",
          "component"
        ),
        g_cols[1]
      )
    )
  )
})

sirius_stats_data <- reactive({
  req(sirius_df(), input$stats_sirius_id_col, input$stats_sirius_col)

  validate(
    need(
      !isTRUE(input$prune_sirius_to_peak) ||
        !is.null(input$file_data),
      paste0(
        "SIRIUS pruning is enabled, but no peak table is uploaded. ",
        "Upload a peak table or switch off SIRIUS pruning."
      )
    )
  )

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

  keep_ids <- selected_peak_ids_for_sirius()

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

selected_peak_ids_for_sirius <- reactive({

  # Switch off means: use all SIRIUS rows.
  if (!isTRUE(input$prune_sirius_to_peak)) {
    return(NULL)
  }

  req(
    raw_df(),
    input$sirius_peak_id_col
  )

  raw <- as.data.frame(
    raw_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  selected_col <- as.character(
    input$sirius_peak_id_col
  )

  ids <- if (
    identical(
      selected_col,
      "__auto__"
    )
  ) {

    as.character(
      seq_len(
        nrow(raw)
      )
    )

  } else {

    validate(
      need(
        selected_col %in% names(raw),
        "Selected peak-table feature ID column was not found."
      )
    )

    trimws(
      as.character(
        raw[[selected_col]]
      )
    )
  }

  ids <- ids[
    !is.na(ids) &
      nzchar(ids)
  ]

  unique(
    ids
  )
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

  out
})

gnps_component_edges <- reactive({

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
    need(
      input$gnps_cluster1_col %in% names(g),
      "GNPS ClusterID1 column was not found."
    ),
    need(
      input$gnps_cluster2_col %in% names(g),
      "GNPS ClusterID2 column was not found."
    ),
    need(
      input$gnps_component_col %in% names(g),
      "GNPS ComponentIndex column was not found."
    )
  )

  edges <- g %>%
    dplyr::transmute(

      ComponentIndex = trimws(
        as.character(
          .data[[input$gnps_component_col]]
        )
      ),

      ClusterID1 = trimws(
        as.character(
          .data[[input$gnps_cluster1_col]]
        )
      ),

      ClusterID2 = trimws(
        as.character(
          .data[[input$gnps_cluster2_col]]
        )
      )
    ) %>%

    dplyr::filter(
      !is.na(ComponentIndex),
      !is.na(ClusterID1),
      !is.na(ClusterID2),
      nzchar(ComponentIndex),
      nzchar(ClusterID1),
      nzchar(ClusterID2),
      ClusterID1 != ClusterID2
    ) %>%

    dplyr::distinct(
      ComponentIndex,
      ClusterID1,
      ClusterID2
    )

  # Optional pruning using the selected peak-table ID column

  edges
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

output$network_component_picker <- renderUI({

  req(
    input$file_gnps_pairs,
    input$stats_selected_class,
    selected_class_component_stats(),
    gnps_component_edges()
  )

  class_stats <- selected_class_component_stats()

  # Remove rows that cannot be displayed as a network
  class_stats <- class_stats %>%
    dplyr::filter(
      !ComponentIndex %in% c(
        "No ComponentIndex match",
        "GNPS pairs not uploaded"
      )
    )

  if (!nrow(class_stats)) {

    return(
      tagList(
        tags$hr(),

        div(
          class = "small-note",
          style = "
            padding: 8px;
            border: 1px solid #dddddd;
            border-radius: 6px;
          ",
          paste0(
            "The selected SIRIUS value does not occur in a ",
            "displayable GNPS component."
          )
        )
      )
    )
  }

  edge_data <- gnps_component_edges()

  # Count total nodes in each available component
  component_totals <- edge_data %>%

    tidyr::pivot_longer(
      cols = c(
        "ClusterID1",
        "ClusterID2"
      ),
      values_to = "ClusterID"
    ) %>%

    dplyr::distinct(
      ComponentIndex,
      ClusterID
    ) %>%

    dplyr::count(
      ComponentIndex,
      name = "Total_nodes"
    )

  component_options <- class_stats %>%

    dplyr::left_join(
      component_totals,
      by = "ComponentIndex"
    ) %>%

    # Exclude components with no surviving edge after pruning
    dplyr::filter(
      !is.na(Total_nodes),
      Total_nodes >= 2
    ) %>%

    dplyr::mutate(
      Display_name = paste0(
        "Component ",
        ComponentIndex,
        " — ",
        Points,
        " selected / ",
        Total_nodes,
        " total"
      )
    ) %>%

    dplyr::arrange(
      dplyr::desc(Points),
      dplyr::desc(Total_nodes),
      ComponentIndex
    )

  if (!nrow(component_options)) {

    return(
      tagList(
        tags$hr(),

        div(
          class = "small-note",
          paste0(
            "No GNPS component containing the selected class ",
"has displayable edges."
          )
        )
      )
    )
  }

  old_component <- isolate(
    input$network_component
  )

  selected_component <- if (
    !is.null(old_component) &&
    old_component %in% component_options$ComponentIndex
  ) {
    old_component
  } else {
    component_options$ComponentIndex[1]
  }

  tagList(

    tags$hr(),

    h3(
      class = "highlight",
      "Interactive component network"
    ),

    pickerInput(
      inputId = "network_component",
      label = "GNPS component to display:",

      choices = stats::setNames(
        component_options$ComponentIndex,
        component_options$Display_name
      ),

      selected = selected_component,
      multiple = FALSE,

      options = list(
        `live-search` = TRUE,
        `style` = "btn-success"
      )
    ),

    radioButtons(
      inputId = "network_label_mode",
      label = "Node labels:",
      choices = c(
        "Selected class only" = "selected",
        "All nodes" = "all",
        "No labels" = "none"
      ),
      selected = "selected"
    ),

    uiOutput(
      "network_volcano_options"
    )
  )
})

output$network_volcano_options <- renderUI({

  if (
    is.null(rv$volcano) ||
    !is.data.frame(rv$volcano) ||
    !nrow(rv$volcano) ||
    !"id" %in% names(rv$volcano)
  ) {

    return(
      div(
        class = "small-note",
        style = "margin-top: 8px;",
        paste0(
          "Run preprocessing to add statistical results from all ",
          "processed comparisons to network-node hover and click details."
        )
      )
    )
  }

  tagList(

    materialSwitch(
      inputId = "network_use_volcano",
      label = "Add all processed statistical comparisons to network nodes",
      value = TRUE,
      status = "success",
      width = "auto"
    ),

    div(
      class = "small-note",
      style = "margin-top: 4px;",
      paste0(
        "For each network node, FC, adjusted p-value, -log10(FDR), ",
        "overall mean, group means, test scale, and significance ",
        "are shown for every processed comparison."
      )
    )
  )
})

selected_component_network_data <- reactive({

  req(
    input$file_gnps_pairs,
    input$network_component,
    input$stats_selected_class,
    gnps_component_edges()
  )

  selected_component <- as.character(
    input$network_component
  )

  component_edges <- gnps_component_edges() %>%

    dplyr::filter(
      ComponentIndex == selected_component
    ) %>%

    dplyr::select(
      ClusterID1,
      ClusterID2
    ) %>%

    dplyr::distinct()

  validate(
    need(
      nrow(component_edges) > 0,
      "The selected component contains no displayable GNPS edges."
    )
  )

  # Build and simplify graph
  graph_object <- igraph::graph_from_data_frame(
    component_edges,
    directed = FALSE
  )

  graph_object <- igraph::simplify(
    graph_object,
    remove.multiple = TRUE,
    remove.loops = TRUE
  )

  node_count <- igraph::vcount(
    graph_object
  )

  validate(
    need(
      node_count <= 500,
      paste0(
        "This component contains ",
        node_count,
        " nodes. Interactive display is limited to 500 nodes."
      )
    )
  )

  # Reproducible force-directed layout
  set.seed(1234)

  coordinates <- igraph::layout_with_fr(
    graph_object,
    niter = 500
  )

  node_names <- igraph::V(
    graph_object
  )$name

  node_data <- tibble::tibble(
    ClusterID = as.character(
      node_names
    ),

    x = as.numeric(
      coordinates[, 1]
    ),

    y = as.numeric(
      coordinates[, 2]
    ),

    Degree = as.integer(
      igraph::degree(
        graph_object
      )
    )
  )

  # Aggregate all values from the selected SIRIUS column
  sirius_node_annotations <- sirius_stats_data() %>%

    dplyr::group_by(
      SIRIUS_ID
    ) %>%

    dplyr::summarise(

      SIRIUS_values = paste(
        sort(
          unique(
            Annotation[
              !is.na(Annotation)
            ]
          )
        ),
        collapse = " | "
      ),

      Selected_class = any(
        Annotation ==
          input$stats_selected_class
      ),

      .groups = "drop"
    ) %>%

    dplyr::rename(
      ClusterID = SIRIUS_ID
    )

  node_data <- node_data %>%

    dplyr::left_join(
      sirius_node_annotations,
      by = "ClusterID"
    ) %>%

    dplyr::mutate(

      Selected_class = dplyr::coalesce(
        Selected_class,
        FALSE
      ),

      Has_SIRIUS = (
        !is.na(SIRIUS_values) &
          nzchar(SIRIUS_values)
      ),

      Node_type = dplyr::case_when(
        Selected_class ~ "Selected class",
        Has_SIRIUS ~ "Other SIRIUS annotation",
        TRUE ~ "No SIRIUS annotation"
      )
    )

  volcano_added <- FALSE

format_network_number <- function(
    x,
    digits = 4
) {

  x <- suppressWarnings(
    as.numeric(x)
  )

  out <- rep(
    "NA",
    length(x)
  )

  ok <- is.finite(x)

  out[ok] <- format(
    signif(
      x[ok],
      digits
    ),
    scientific = FALSE,
    trim = TRUE
  )

  out
}


# Empty structure in case volcano statistics are unavailable
volcano_stats_long <- tibble::tibble(

  ClusterID = character(),
  Comparison = character(),

  Group_num = character(),
  Group_den = character(),

  FC = numeric(),
  Adj_p = numeric(),
  Adj_p_log = numeric(),

  Mean = numeric(),
  mean_num = numeric(),
  mean_den = numeric(),

  TestScale = character(),
  Significant_default = character(),

  GNPS_annotation = character(),

  Comparison_hover = character()
)


if (
  isTRUE(input$network_use_volcano) &&
  !is.null(rv$volcano) &&
  is.data.frame(rv$volcano) &&
  nrow(rv$volcano) &&
  "id" %in% names(rv$volcano)
) {

  volc <- as.data.frame(
    rv$volcano,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


  n_volc <- nrow(volc)


  get_chr <- function(column_name) {

    if (!column_name %in% names(volc)) {

      return(
        rep(
          NA_character_,
          n_volc
        )
      )
    }

    as.character(
      volc[[column_name]]
    )
  }


  get_num <- function(column_name) {

    if (!column_name %in% names(volc)) {

      return(
        rep(
          NA_real_,
          n_volc
        )
      )
    }

    suppressWarnings(
      as.numeric(
        volc[[column_name]]
      )
    )
  }


  comparison_values <- get_chr(
    "Groups"
  )

  comparison_values[
    is.na(comparison_values) |
      !nzchar(
        trimws(
          comparison_values
        )
      )
  ] <- "Comparison"


  gnps_values <- if (
    "GNPS_annotation" %in% names(volc)
  ) {

    clean_missing_text(
      volc$GNPS_annotation
    )

  } else {

    rep(
      NA_character_,
      n_volc
    )
  }


  volcano_stats_long <- tibble::tibble(

    ClusterID = trimws(
      as.character(
        volc$id
      )
    ),

    Comparison = comparison_values,

    Group_num = get_chr(
      "Group_num"
    ),

    Group_den = get_chr(
      "Group_den"
    ),

    FC = get_num(
      "FC"
    ),

    Adj_p = get_num(
      "Adj.p-value"
    ),

    Adj_p_log = get_num(
      "Adj.p-value.log"
    ),

    Mean = get_num(
      "Mean"
    ),

    mean_num = get_num(
      "mean_num"
    ),

    mean_den = get_num(
      "mean_den"
    ),

    TestScale = get_chr(
      "TestScale"
    ),

    Significant_default = get_chr(
      "Significant_default"
    ),

    GNPS_annotation = gnps_values
  ) %>%

    dplyr::filter(
      !is.na(ClusterID),
      nzchar(ClusterID)
    ) %>%

    dplyr::distinct(
      ClusterID,
      Comparison,
      .keep_all = TRUE
    )


  if (nrow(volcano_stats_long)) {

    group_num_label <- volcano_stats_long$Group_num
    group_den_label <- volcano_stats_long$Group_den


    group_num_label[
      is.na(group_num_label) |
        !nzchar(group_num_label)
    ] <- "Group 1"


    group_den_label[
      is.na(group_den_label) |
        !nzchar(group_den_label)
    ] <- "Group 2"


    test_scale_label <- volcano_stats_long$TestScale

    test_scale_label[
      is.na(test_scale_label) |
        !nzchar(test_scale_label)
    ] <- "NA"


    significant_label <-
      volcano_stats_long$Significant_default

    significant_label[
      is.na(significant_label) |
        !nzchar(significant_label)
    ] <- "NA"


    volcano_stats_long$Comparison_hover <- paste0(

      "<br><br><b>Comparison: ",
      htmltools::htmlEscape(
        volcano_stats_long$Comparison
      ),
      "</b>",

      "<br>FC: ",
      format_network_number(
        volcano_stats_long$FC
      ),

      "<br>Adjusted p-value: ",
      format_network_number(
        volcano_stats_long$Adj_p
      ),

      "<br>-log10(FDR): ",
      format_network_number(
        volcano_stats_long$Adj_p_log
      ),

      "<br>Mean intensity: ",
      format_network_number(
        volcano_stats_long$Mean
      ),

      "<br>Mean [",
      htmltools::htmlEscape(
        group_num_label
      ),
      "]: ",
      format_network_number(
        volcano_stats_long$mean_num
      ),

      "<br>Mean [",
      htmltools::htmlEscape(
        group_den_label
      ),
      "]: ",
      format_network_number(
        volcano_stats_long$mean_den
      ),

      "<br>Test scale: ",
      htmltools::htmlEscape(
        test_scale_label
      ),

      "<br>Significant: ",
      htmltools::htmlEscape(
        significant_label
      )
    )


    statistics_by_node <- volcano_stats_long %>%

      dplyr::group_by(
        ClusterID
      ) %>%

      dplyr::summarise(

        Statistical_results = paste0(
          Comparison_hover,
          collapse = ""
        ),

        .groups = "drop"
      )


    gnps_by_node <- volcano_stats_long %>%

      dplyr::group_by(
        ClusterID
      ) %>%

      dplyr::summarise(

        GNPS_annotation = {

          values <- clean_missing_text(
            GNPS_annotation
          )

          values <- unique(
            values[
              !is.na(values)
            ]
          )

          if (length(values)) {

            paste(
              values,
              collapse = " | "
            )

          } else {

            NA_character_
          }
        },

        .groups = "drop"
      )


    node_data <- node_data %>%

      dplyr::left_join(
        statistics_by_node,
        by = "ClusterID"
      ) %>%

      dplyr::left_join(
        gnps_by_node,
        by = "ClusterID"
      )


    volcano_added <- TRUE
  }
}


if (!"Statistical_results" %in% names(node_data)) {
  node_data$Statistical_results <- NA_character_
}


if (!"GNPS_annotation" %in% names(node_data)) {
  node_data$GNPS_annotation <- NA_character_
}


sirius_text <- node_data$SIRIUS_values

sirius_text[
  is.na(sirius_text) |
    !nzchar(sirius_text)
] <- "NA"


gnps_text <- node_data$GNPS_annotation

gnps_text[
  is.na(gnps_text) |
    !nzchar(gnps_text)
] <- "NA"


node_data$Hover <- paste0(

  "<b>Cluster ID:</b> ",
  node_data$ClusterID,

  "<br><b>ComponentIndex:</b> ",
  selected_component,

  "<br><b>Node degree:</b> ",
  node_data$Degree,

  "<br><b>",
  input$stats_sirius_col,
  ":</b> ",
  sirius_text,

  "<br><b>Selected class:</b> ",
  ifelse(
    node_data$Selected_class,
    "Yes",
    "No"
  )
)


if (volcano_added) {

  stats_hover <-
    node_data$Statistical_results

  stats_hover[
    is.na(stats_hover) |
      !nzchar(stats_hover)
  ] <- paste0(
    "<br><br>",
    "No processed statistical results matched this node."
  )


  node_data$Hover <- paste0(

    node_data$Hover,

    "<br><b>GNPS annotation:</b> ",
    gnps_text,

    "<br><br><b>Processed statistical comparisons</b>",

    stats_hover
  )
}

  label_mode <- input$network_label_mode %||%
    "selected"

  node_data$Node_label <- dplyr::case_when(
    label_mode == "all" ~ node_data$ClusterID,
    label_mode == "selected" &
      node_data$Selected_class ~ node_data$ClusterID,
    TRUE ~ ""
  )

  # Create edge coordinates using the simplified graph
  simplified_edges <- igraph::as_data_frame(
    graph_object,
    what = "edges"
  ) %>%

    dplyr::transmute(
      ClusterID1 = as.character(from),
      ClusterID2 = as.character(to)
    )

  edge_coordinates <- simplified_edges %>%

    dplyr::left_join(
      node_data %>%
        dplyr::select(
          ClusterID,
          x1 = x,
          y1 = y
        ),
      by = c(
        "ClusterID1" = "ClusterID"
      )
    ) %>%

    dplyr::left_join(
      node_data %>%
        dplyr::select(
          ClusterID,
          x2 = x,
          y2 = y
        ),
      by = c(
        "ClusterID2" = "ClusterID"
      )
    )

  list(
  component = selected_component,
  graph = graph_object,
  nodes = node_data,
  edges = edge_coordinates,
  volcano_stats = volcano_stats_long,
  volcano_added = volcano_added,
  node_count = igraph::vcount(graph_object),
  edge_count = igraph::ecount(graph_object)
)
})

output$gnps_component_network <- plotly::renderPlotly({

  network <- selected_component_network_data()

  node_data <- network$nodes
  edge_data <- network$edges

  plot_object <- plotly::plot_ly(
    source = "gnps_component_network"
  )

  # GNPS edges
  plot_object <- plot_object %>%

    plotly::add_segments(
      data = edge_data,

      x = ~x1,
      y = ~y1,
      xend = ~x2,
      yend = ~y2,

      inherit = FALSE,

      line = list(
        color = "rgba(130,130,130,0.55)",
        width = 1
      ),

      hoverinfo = "skip",
      showlegend = FALSE
    )

  node_types <- c(
    "Selected class",
    "Other SIRIUS annotation",
    "No SIRIUS annotation"
  )

  node_colors <- c(
    "Selected class" = "#E74C3C",
    "Other SIRIUS annotation" = "#66CDAA",
    "No SIRIUS annotation" = "#BDBDBD"
  )

  node_sizes <- c(
    "Selected class" = 16,
    "Other SIRIUS annotation" = 11,
    "No SIRIUS annotation" = 9
  )

  for (current_type in node_types) {

    current_nodes <- node_data %>%
      dplyr::filter(
        Node_type == current_type
      )

    if (!nrow(current_nodes)) {
      next
    }

    plot_object <- plot_object %>%

      plotly::add_trace(
        data = current_nodes,

        x = ~x,
        y = ~y,

        type = "scatter",
        mode = "markers+text",

        text = ~Node_label,
        hovertext = ~Hover,
        hoverinfo = "text",

        key = ~ClusterID,

        name = current_type,
        inherit = FALSE,

        textposition = "top center",

        textfont = list(
          size = 10,
          color = "#222222"
        ),

        marker = list(
          size = unname(
            node_sizes[
              current_type
            ]
          ),

          color = unname(
            node_colors[
              current_type
            ]
          ),

          line = list(
            color = "#333333",
            width = 1
          )
        )
      )
  }

  plot_object %>%

    plotly::layout(

      title = list(
        text = paste0(
          "GNPS Component ",
          network$component,
          " — ",
          network$node_count,
          " nodes, ",
          network$edge_count,
          " edges"
        )
      ),

      xaxis = list(
        visible = FALSE,
        showgrid = FALSE,
        zeroline = FALSE
      ),

      yaxis = list(
        visible = FALSE,
        showgrid = FALSE,
        zeroline = FALSE,
        scaleanchor = "x",
        scaleratio = 1
      ),

      hovermode = "closest",
      dragmode = "pan",

      legend = list(
        orientation = "h",
        x = 0,
        y = -0.05
      ),

      margin = list(
        l = 20,
        r = 20,
        b = 70,
        t = 70
      )
    ) %>%

    plotly::config(
      displaylogo = FALSE,

      modeBarButtonsToRemove = c(
        "lasso2d",
        "select2d"
      )
    )
})

output$gnps_network_node_details <- renderUI({

  click <- plotly::event_data(
    "plotly_click",
    source = "gnps_component_network"
  )


  if (
    is.null(click) ||
    is.null(click$key) ||
    !length(click$key)
  ) {

    return(
      div(
        class = "small-note",
        style = "margin-top: 8px;",
        "Click a network node to show its information."
      )
    )
  }


  network <- selected_component_network_data()


  clicked_id <- as.character(
    click$key[1]
  )


  row <- network$nodes %>%

    dplyr::filter(
      ClusterID == clicked_id
    ) %>%

    dplyr::slice(1)


  if (!nrow(row)) {
    return(NULL)
  }


  sirius_value <- row$SIRIUS_values[1]


  if (
    is.na(sirius_value) ||
    !nzchar(sirius_value)
  ) {

    sirius_value <- "NA"
  }


  gnps_value <- row$GNPS_annotation[1]


  if (
    is.na(gnps_value) ||
    !nzchar(gnps_value)
  ) {

    gnps_value <- "NA"
  }


  format_one <- function(x) {

    x <- suppressWarnings(
      as.numeric(x)
    )

    if (
      !length(x) ||
      !is.finite(x[1])
    ) {
      return("NA")
    }

    format(
      signif(
        x[1],
        5
      ),
      scientific = FALSE,
      trim = TRUE
    )
  }


  stats_rows <- network$volcano_stats %>%

    dplyr::filter(
      ClusterID == clicked_id
    )


  stats_table_ui <- if (!nrow(stats_rows)) {

    div(
      class = "small-note",
      style = "margin-top: 10px;",
      "No processed statistical comparison matched this network node."
    )

  } else {


    table_rows <- lapply(

      seq_len(
        nrow(stats_rows)
      ),

      function(i) {


        current <- stats_rows[
          i,
          ,
          drop = FALSE
        ]


        group_num <- current$Group_num[1]
        group_den <- current$Group_den[1]


        if (
          is.na(group_num) ||
          !nzchar(group_num)
        ) {
          group_num <- "Group 1"
        }


        if (
          is.na(group_den) ||
          !nzchar(group_den)
        ) {
          group_den <- "Group 2"
        }


        test_scale <- current$TestScale[1]

        if (
          is.na(test_scale) ||
          !nzchar(test_scale)
        ) {
          test_scale <- "NA"
        }


        significant_value <-
          current$Significant_default[1]

        if (
          is.na(significant_value) ||
          !nzchar(significant_value)
        ) {
          significant_value <- "NA"
        }


        tags$tr(

          tags$td(
            current$Comparison[1]
          ),

          tags$td(
            group_num
          ),

          tags$td(
            format_one(
              current$mean_num
            )
          ),

          tags$td(
            group_den
          ),

          tags$td(
            format_one(
              current$mean_den
            )
          ),

          tags$td(
            format_one(
              current$FC
            )
          ),

          tags$td(
            format_one(
              current$Adj_p
            )
          ),

          tags$td(
            format_one(
              current$Adj_p_log
            )
          ),

          tags$td(
            format_one(
              current$Mean
            )
          ),

          tags$td(
            test_scale
          ),

          tags$td(
            significant_value
          )
        )
      }
    )


    tagList(

      h4(
        style = "margin-top: 14px;",
        "All processed statistical comparisons"
      ),


      div(
        style = "overflow-x:auto;",


        tags$table(

          class =
            "table table-striped table-bordered table-condensed",


          tags$thead(

            tags$tr(

              tags$th("Comparison"),

              tags$th("Group 1"),

              tags$th("Mean 1"),

              tags$th("Group 2"),

              tags$th("Mean 2"),

              tags$th("FC"),

              tags$th("Adjusted p-value"),

              tags$th("-log10(FDR)"),

              tags$th("Mean intensity"),

              tags$th("Test scale"),

              tags$th("Significant")
            )
          ),


          tags$tbody(
            table_rows
          )
        )
      )
    )
  }


  div(

    style = "
      background:#ffffffcc;
      padding:10px;
      border-radius:10px;
      border:1px solid #dddddd;
      margin-top:10px;
    ",


    h4(
      paste0(
        "Selected network node: ",
        clicked_id
      )
    ),


    tags$ul(

      tags$li(
        strong("ComponentIndex: "),
        network$component
      ),

      tags$li(
        strong("Node degree: "),
        row$Degree[1]
      ),

      tags$li(
        strong(
          paste0(
            input$stats_sirius_col,
            ": "
          )
        ),
        sirius_value
      ),

      tags$li(
        strong("Matches selected class: "),
        ifelse(
          isTRUE(
            row$Selected_class[1]
          ),
          "Yes",
          "No"
        )
      ),

      tags$li(
        strong("GNPS annotation: "),
        gnps_value
      )
    ),


    stats_table_ui
  )
})

output$gnps_network_section <- renderUI({

  if (is.null(input$file_gnps_pairs)) {
    return(NULL)
  }

  if (
    is.null(input$network_component) ||
    !nzchar(
      as.character(
        input$network_component
      )
    )
  ) {

    return(
      tagList(
        tags$hr(),

        div(
          class = "small-note",
          paste0(
            "Choose a SIRIUS class that occurs in a GNPS ",
            "component to display its network."
          )
        )
      )
    )
  }

  tagList(

    tags$hr(),

    h3(
      "Interactive selected GNPS component"
    ),

    div(
      class = "small-note",
      style = "margin-bottom: 8px;",

      HTML(
        paste0(
          "<b>Red:</b> selected SIRIUS class &nbsp;|&nbsp; ",
          "<b>Green:</b> another SIRIUS value &nbsp;|&nbsp; ",
          "<b>Grey:</b> no value in the selected SIRIUS column"
        )
      )
    ),

    withSpinner(
      plotlyOutput(
        "gnps_component_network",
        height = "650px"
      ),
      type = 8,
      color = "#66CDAA"
    ),

    uiOutput(
      "gnps_network_node_details"
    )
  )
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
    withSpinner(DTOutput("selected_class_component_table"), type = 8, color = "#66CDAA"),
    uiOutput("gnps_network_section")
  )
})

  # ---- Process button ----
  observeEvent(input$run_proc, {
  req(built(), labels_vec())
  labs_pre <- labels_vec()
  if (stop_if_one_group(labs_pre)) {
    return(NULL)
  }
  comparison_mode <- input$comparison_mode %||% "reference"
  manual_pairs <- NULL
  if (identical(comparison_mode, "manual")) {
    if (
      is.null(input$manual_comparisons) ||
      length(input$manual_comparisons) == 0
    ) {
      showNotification(
        "Select at least one manual comparison.",
        type = "error",
        duration = 6
      )
      return(NULL)
    }
    manual_pairs <- selected_manual_comparisons()
  } else {
    req(input$ref_group)
  }
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
        adj = input$p_adjust %||% "BH",
        paired = isTRUE(input$paired),
        eqvar = isTRUE(input$eqvar),
        pseudocount = 1.1,
        log2_test = isTRUE(input$log2_test),
        scale_data = isTRUE(input$standard_scaling),

        ref_group = if (
          identical(comparison_mode, "reference")
        ) {
          input$ref_group
        } else {
          NULL
        },

        comparisons = manual_pairs
      )

      incProgress(0.05, detail = "Joining mz/rt/id")
      volc <- volc %>% left_join(fmap, by = "Feature")

      incProgress(
  0.02,
  detail = "Joining additional peak-table columns"
)

if (isTRUE(input$use_peak_extra_cols)) {

  selected_peak_cols <- intersect(
    input$peak_extra_cols %||% character(0),
    names(built()$raw)
  )

  if (length(selected_peak_cols)) {

    peak_raw <- as.data.frame(
      built()$raw,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    validate(
      need(
        nrow(peak_raw) == nrow(fmap),
        paste0(
          "Peak-table metadata could not be joined because ",
          "the peak table and feature map have different row counts."
        )
      )
    )

    peak_colmap <- make_prefixed_colmap(
      selected_peak_cols,
      prefix = "Peak_"
    )

    peak_extra <- peak_raw[
      ,
      names(peak_colmap),
      drop = FALSE
    ]

    names(peak_extra) <- unname(
      peak_colmap
    )

    peak_extra$Feature <- as.character(
      fmap$Feature
    )

    peak_extra <- peak_extra[
      ,
      c(
        "Feature",
        unname(peak_colmap)
      ),
      drop = FALSE
    ]

    volc <- volc %>%
      dplyr::left_join(
        peak_extra,
        by = "Feature"
      )
  }
}

      # default annotation columns
      volc$`NPC#class` <- NA_character_
      volc$`ClassyFire#class` <- NA_character_
      volc$GNPS_annotation <- NA_character_

      incProgress(
      0.05,
      detail = "Joining SIRIUS (optional)"
    )

if (
  isTRUE(input$use_sirius) &&
  !is.null(input$file_sirius)
) {

  s <- as.data.frame(
    sirius_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  req(
    input$sirius_idcol,
    input$sirius_npcol,
    input$sirius_cfcol
  )

  validate(
    need(
      input$sirius_idcol %in% names(s),
      "Selected SIRIUS Feature ID column was not found."
    ),

    need(
      input$sirius_npcol %in% names(s),
      "Selected SIRIUS NPC column was not found."
    ),

    need(
      input$sirius_cfcol %in% names(s),
      "Selected SIRIUS ClassyFire column was not found."
    )
  )

  selected_sirius_extra <- character(0)

  if (isTRUE(input$use_sirius_extra_cols)) {

    selected_sirius_extra <- intersect(
      input$sirius_extra_cols %||% character(0),
      names(s)
    )

    selected_sirius_extra <- setdiff(
      selected_sirius_extra,
      c(
        input$sirius_idcol,
        input$sirius_npcol,
        input$sirius_cfcol
      )
    )
  }

  sirius_colmap <- make_prefixed_colmap(
    selected_sirius_extra,
    prefix = "SIRIUS_"
  )

  # Build the core SIRIUS join table
  ss <- data.frame(
    id = trimws(
      as.character(
        s[[input$sirius_idcol]]
      )
    ),

    `NPC#class` = clean_missing_text(
      s[[input$sirius_npcol]]
    ),

    `ClassyFire#class` = clean_missing_text(
      s[[input$sirius_cfcol]]
    ),

    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Add selected additional SIRIUS columns.
  # Numeric columns remain numeric.
  if (length(sirius_colmap)) {

    for (original_name in names(sirius_colmap)) {

      output_name <- sirius_colmap[[
        original_name
      ]]

      value <- s[[
        original_name
      ]]

      if (
        is.character(value) ||
        is.factor(value)
      ) {
        value <- clean_missing_text(
          value
        )
      }

      ss[[
        output_name
      ]] <- value
    }
  }

  ss <- ss %>%
    dplyr::filter(
      !is.na(id),
      nzchar(id)
    )

  # Prevent duplicate SIRIUS IDs from multiplying volcano rows.
  # The first row is retained if duplicate IDs are present.
  duplicate_count <- sum(
    duplicated(ss$id)
  )

  if (duplicate_count > 0) {

    showNotification(
      paste0(
        "SIRIUS contains ",
        duplicate_count,
        " duplicated mapping ID row(s). ",
        "The first row for each ID was retained."
      ),
      type = "warning",
      duration = 7
    )

    ss <- ss %>%
      dplyr::distinct(
        id,
        .keep_all = TRUE
      )
  }

peak_table_ids <- unique(
  trimws(
    as.character(
      volc$id
    )
  )
)

sirius_ids <- unique(
  trimws(
    as.character(
      ss$id
    )
  )
)

peak_table_ids <- peak_table_ids[
  !is.na(peak_table_ids) &
    nzchar(peak_table_ids)
]

sirius_ids <- sirius_ids[
  !is.na(sirius_ids) &
    nzchar(sirius_ids)
]

matched_sirius_ids <- sum(
  peak_table_ids %in% sirius_ids
)

volc <- volc %>%

  dplyr::mutate(
    id = trimws(
      as.character(id)
    )
  ) %>%

  dplyr::left_join(
    ss,
    by = "id",
    suffix = c("", ".sirius")
  ) %>%

  dplyr::mutate(

    `NPC#class` = dplyr::coalesce(
      .data[["NPC#class.sirius"]],
      .data[["NPC#class"]]
    ),

    `ClassyFire#class` = dplyr::coalesce(
      .data[["ClassyFire#class.sirius"]],
      .data[["ClassyFire#class"]]
    )
  ) %>%

  dplyr::select(
    -dplyr::any_of(
      c(
        "NPC#class.sirius",
        "ClassyFire#class.sirius"
      )
    )
  )

if (matched_sirius_ids == 0) {

  showNotification(
    paste0(
      "SIRIUS summary was uploaded, but no IDs matched. ",
      "Check that the peak-table Row ID column and selected ",
      "SIRIUS Feature ID column contain the same values."
    ),
    type = "warning",
    duration = 8
  )

} else {

  showNotification(
    paste0(
      "SIRIUS annotation joined successfully: ",
      format(
        matched_sirius_ids,
        big.mark = ",",
        scientific = FALSE
      ),
      " of ",
      format(
        length(peak_table_ids),
        big.mark = ",",
        scientific = FALSE
      ),
      " unique peak-table IDs matched."
    ),
    type = "message",
    duration = 6
  )
}
}

incProgress(
  0.05,
  detail = "Joining GNPS annotation (optional)"
)

if (
  isTRUE(input$use_gnps_annotation) &&
  !is.null(input$file_gnps_annotation)
) {

  g <- as.data.frame(
    gnps_annotation_df(),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  req(
    input$gnps_annotation_idcol,
    input$gnps_annotation_col
  )

  validate(
    need(
      input$gnps_annotation_idcol %in% names(g),
      "Selected GNPS ID column was not found."
    ),

    need(
      input$gnps_annotation_col %in% names(g),
      "Selected primary GNPS annotation column was not found."
    )
  )

  selected_gnps_extra <- character(0)

  if (isTRUE(input$use_gnps_extra_cols)) {

    selected_gnps_extra <- intersect(
      input$gnps_extra_cols %||% character(0),
      names(g)
    )

    selected_gnps_extra <- setdiff(
      selected_gnps_extra,
      c(
        input$gnps_annotation_idcol,
        input$gnps_annotation_col
      )
    )
  }

  gnps_colmap <- make_prefixed_colmap(
    selected_gnps_extra,
    prefix = "GNPS_"
  )

  gnps_primary <- data.frame(
    id = trimws(
      as.character(
        g[[input$gnps_annotation_idcol]]
      )
    ),

    GNPS_annotation = clean_missing_text(
      g[[input$gnps_annotation_col]]
    ),

    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%

    dplyr::filter(
      !is.na(id),
      nzchar(id)
    ) %>%

    dplyr::group_by(id) %>%

    dplyr::summarise(

      GNPS_annotation = {

        values <- unique(
          GNPS_annotation[
            !is.na(GNPS_annotation)
          ]
        )

        if (length(values)) {
          paste(
            values,
            collapse = " | "
          )
        } else {
          NA_character_
        }
      },

      .groups = "drop"
    )

  gnps_extra <- NULL

  if (length(gnps_colmap)) {

    gnps_extra <- data.frame(
      id = trimws(
        as.character(
          g[[input$gnps_annotation_idcol]]
        )
      ),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    for (original_name in names(gnps_colmap)) {

      output_name <- gnps_colmap[[
        original_name
      ]]

      value <- g[[
        original_name
      ]]

      # Convert empty character values to real NA
      if (
        is.character(value) ||
        is.factor(value)
      ) {
        value <- clean_missing_text(
          value
        )
      }

      # Numeric GNPS columns remain numeric
      gnps_extra[[
        output_name
      ]] <- value
    }

    gnps_extra <- gnps_extra %>%
      dplyr::filter(
        !is.na(id),
        nzchar(id)
      )

    duplicate_extra_ids <- sum(
      duplicated(gnps_extra$id)
    )

    if (duplicate_extra_ids > 0) {

      showNotification(
        paste0(
          "GNPS contains ",
          duplicate_extra_ids,
          " duplicated ID row(s). ",
          "The primary annotations were combined, while the first row ",
          "was retained for each additional GNPS column."
        ),
        type = "warning",
        duration = 8
      )
    }

    gnps_extra <- gnps_extra %>%
      dplyr::distinct(
        id,
        .keep_all = TRUE
      )
  }

  # Combine primary and additional GNPS columns
  gnps_join <- gnps_primary

  if (!is.null(gnps_extra)) {

    gnps_join <- gnps_join %>%
      dplyr::left_join(
        gnps_extra,
        by = "id"
      )
  }

  peak_table_ids <- unique(
  trimws(
    as.character(
      volc$id
    )
  )
)

peak_table_ids <- peak_table_ids[
  !is.na(peak_table_ids) &
    nzchar(peak_table_ids)
]

gnps_ids <- unique(
  trimws(
    as.character(
      gnps_join$id
    )
  )
)

gnps_ids <- gnps_ids[
  !is.na(gnps_ids) &
    nzchar(gnps_ids)
]

matched_gnps_ids <- sum(
  peak_table_ids %in% gnps_ids
)

  volc <- volc %>%

    dplyr::mutate(
      id = trimws(
        as.character(id)
      )
    ) %>%

    dplyr::left_join(
      gnps_join,
      by = "id",
      suffix = c("", ".gnps")
    ) %>%

    dplyr::mutate(
      GNPS_annotation = dplyr::coalesce(
        .data[["GNPS_annotation.gnps"]],
        .data[["GNPS_annotation"]]
      )
    ) %>%

    dplyr::select(
      -dplyr::any_of(
        "GNPS_annotation.gnps"
      )
    )

  if (matched_gnps_ids == 0) {

    showNotification(
      paste0(
        "GNPS annotation was uploaded, but no IDs matched. ",
        "Check that the peak-table Row ID column and selected ",
        "GNPS ID column contain the same values."
      ),
      type = "warning",
      duration = 8
    )

  } else {

    showNotification(
  paste0(
    "GNPS annotation joined successfully: ",
    format(
      matched_gnps_ids,
      big.mark = ",",
      scientific = FALSE
    ),
    " of ",
    format(
      length(peak_table_ids),
      big.mark = ",",
      scientific = FALSE
    ),
    " unique peak-table IDs matched."
  ),
  type = "message",
  duration = 6
)
  }
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
      updateRadioButtons(session, "volcano_y_axis", selected = "fdr")
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

  x <- clean_missing_text(
    x
  )

  x <- x[
    !is.na(x) &
      nzchar(x)
  ]

  sort(
    unique(x)
  )
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
        "volcano_y_axis",
        "Y-axis:",
        choices = c(
          "-log10(FDR)" = "fdr",
          "Mean log10(Intensity)" = "mean"
        ),
        selected = "fdr",
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
      uiOutput("selected_feature_panel")
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
    if (!is.finite(pcut) || pcut < 0 || pcut > 1) pcut <- 0.05

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

    use_mean_y <- identical(
        input$volcano_y_axis %||% "fdr",
        "mean"
      )

      dd$plot_y <- if (use_mean_y) {
        log10(dd$Mean + 1.1)
      } else {
        dd$`Adj.p-value.log`
      }

      y_title <- if (use_mean_y) {
        "Mean log10(Intensity)"
      } else {
        "-log10(FDR)"
      }

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
      "<br>ClassyFire: ", dd$`ClassyFire#class`,
      "<br>GNPS annotation: ", dd$GNPS_annotation
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
  )
)

# Add the FDR cutoff line only when FDR is used as y-axis
if (!use_mean_y) {
  shapes <- c(
    shapes,
    list(
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
  )
}

    if (input$color_by == "Groups") {
      plot_ly(
        data = dd,
        x = ~FC, y = ~plot_y,
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
          yaxis = list(title = y_title),
          legend = list(title = list(text = "Comparison"))
        ) %>% event_register("plotly_click")
    } else {
      plot_ly(
        data = dd,
        x = ~FC, y = ~plot_y,
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
          yaxis = list(title = y_title),
          legend = list(title = list(text = "Comparison"))
        ) %>% event_register("plotly_click")
    }
  })

  # ---- Feature plot on click
  output$selected_feature_panel <- renderUI({
  req(procReady(), rv$volcano)

  click <- plotly::event_data(
    "plotly_click",
    source = "volcano"
  )

  if (
    is.null(click) ||
    is.null(click$key) ||
    !length(click$key)
  ) {
    return(NULL)
  }

  key <- as.character(click$key[[1]])

  parts <- stringr::str_split_fixed(
    key,
    "__",
    2
  )

  comp <- parts[1, 1]
  feat <- parts[1, 2]

  row <- rv$volcano %>%
    dplyr::filter(
      Groups == comp,
      Feature == feat
    ) %>%
    dplyr::slice(1)

  if (nrow(row) == 0) {
    return(NULL)
  }

  # Prefer the original peak-table ID.
  # Fall back to the internal Feature value.
  feature_text <- feat

  if (
    "id" %in% names(row) &&
    !is.na(row$id[[1]]) &&
    nzchar(trimws(as.character(row$id[[1]])))
  ) {
    feature_text <- as.character(row$id[[1]])
  }

  mz_value <- suppressWarnings(
    as.numeric(row$mz[[1]])
  )

  rt_value <- suppressWarnings(
    as.numeric(row$RT[[1]])
  )

  mz_text <- if (is.finite(mz_value)) {
    format(
      mz_value,
      digits = 12,
      scientific = FALSE,
      trim = TRUE
    )
  } else {
    "NA"
  }

  rt_text <- if (is.finite(rt_value)) {
    format(
      rt_value,
      digits = 8,
      scientific = FALSE,
      trim = TRUE
    )
  } else {
    "NA"
  }

  peak_detail_map <- make_prefixed_colmap(

  if (isTRUE(input$use_peak_extra_cols)) {
    input$peak_extra_cols %||% character(0)
  } else {
    character(0)
  },

  prefix = "Peak_"
)

sirius_detail_map <- make_prefixed_colmap(

  if (
    isTRUE(input$use_sirius) &&
    isTRUE(input$use_sirius_extra_cols)
  ) {
    input$sirius_extra_cols %||% character(0)
  } else {
    character(0)
  },

  prefix = "SIRIUS_"
)

gnps_detail_map <- make_prefixed_colmap(

  if (
    isTRUE(input$use_gnps_annotation) &&
    isTRUE(input$use_gnps_extra_cols)
  ) {
    input$gnps_extra_cols %||% character(0)
  } else {
    character(0)
  },

  prefix = "GNPS_"
)

detail_labels <- c(

  stats::setNames(
    paste0(
      "Peak table: ",
      names(peak_detail_map)
    ),
    unname(peak_detail_map)
  ),

  stats::setNames(
    paste0(
      "SIRIUS: ",
      names(sirius_detail_map)
    ),
    unname(sirius_detail_map)
  ),

  stats::setNames(
    paste0(
      "GNPS: ",
      names(gnps_detail_map)
    ),
    unname(gnps_detail_map)
  )
)

detail_cols <- intersect(
  names(detail_labels),
  names(row)
)

additional_details_ui <- NULL

if (length(detail_cols)) {

  additional_details_ui <- div(

    tags$hr(),

    tags$h5(
      style = "
        font-weight: bold;
        color: #2c3e50;
        margin-bottom: 10px;
      ",
      "Additional information"
    ),

    lapply(
      detail_cols,
      function(column_name) {

        div(
          style = "
            display: grid;
            grid-template-columns: minmax(180px, 35%) 1fr;
            gap: 10px;
            padding: 5px 0;
            border-bottom: 1px solid #eeeeee;
            overflow-wrap: anywhere;
          ",

          tags$strong(
            paste0(
              detail_labels[[column_name]],
              ":"
            )
          ),

          tags$span(
            format_extra_value(
              row[[column_name]]
            )
          )
        )
      }
    )
  )
}

  tagList(
    div(
      style = "
        background-color: white;
        border: 2px solid #66CDAA;
        border-radius: 8px;
        padding: 12px 15px;
        margin-bottom: 12px;
      ",

      tags$h4(
        style = "
          margin-top: 0;
          margin-bottom: 12px;
          font-weight: bold;
          color: #2c3e50;
        ",
        "Feature ID: ",
        tags$span(
          style = "color: #18bc9c;",
          feature_text
        )
      ),

      fluidRow(
        column(
          width = 6,

          tags$label(
            `for` = "selected_mz_text",
            style = "font-weight: bold;",
            "m/z:"
          ),

          tags$input(
            id = "selected_mz_text",
            type = "text",
            class = "form-control",
            value = mz_text,
            readonly = "readonly",

            # Clicking or focusing selects the complete value
            onclick = "this.select();",
            onfocus = "this.select();"
          )
        ),

        column(
          width = 6,

          tags$label(
            `for` = "selected_rt_text",
            style = "font-weight: bold;",
            "RT:"
          ),

          tags$input(
            id = "selected_rt_text",
            type = "text",
            class = "form-control",
            value = rt_text,
            readonly = "readonly",

            onclick = "this.select();",
            onfocus = "this.select();"
          )
        )
      ),

      additional_details_ui
    ),

    plotlyOutput(
      "feature_plot",
      height = "260px"
    )
  )
})

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
