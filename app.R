#Sys.setlocale("LC_ALL", "English_United States.1252")
#Sys.setenv(LANG = "en_US.UTF-8")
#Sys.setlocale("LC_ALL", "C.UTF-8")   # Use a universally available UTF-8 locale
#options(encoding = "UTF-8")
#rsconnect::deployApp()

# app.R --------------------------------------------------------------------
# Volcano Explorer (2 pages):
#   1) Load & Process
#   2) Volcano explorer
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinyWidgets)
  library(shinythemes)
  library(DT)
  library(vroom)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(stringr)
  library(plotly)
  library(MsCoreUtils)   # bpca imputation
  library(viridisLite)   # colors for Mean mode
  library(tibble)
  library(pcaMethods)
})

# ----------------------------- Helpers -----------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

clean_mzmine_export <- function(df) {
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) > 0) {
    last <- df[[ncol(df)]]
    if (all(is.na(last)) || all(trimws(as.character(last)) == "")) {
      if (grepl("^Unnamed", names(df)[ncol(df)])) {
        df <- df[, -ncol(df), drop = FALSE]
      }
    }
  }
  df
}

multi_sample_idx <- function(cols, kws) {
  kws <- as.character(kws)
  kws <- kws[nzchar(kws)]
  if (!length(kws)) return(integer(0))
  hits <- Reduce(`|`, lapply(kws, function(k) grepl(k, cols, fixed = TRUE)))
  which(hits)
}

clean_sample_names <- function(x) {
  x <- gsub(" Peak area$", "", x, ignore.case = TRUE)
  x <- gsub("\\.(mzML|mzXML|raw)$", "", x, ignore.case = TRUE)
  x <- trimws(x)
  x
}

stop_if_one_group <- function(labs) {
  labs <- trimws(as.character(labs))
  labs <- labs[nzchar(labs)]
  u <- unique(labs)

  if (length(u) < 2) {
    shiny::showNotification(
      paste0("Need at least 2 different groups in 'Label' to run statistics."),
      type = "error",
      duration = 8
    )
    return(TRUE)  
  }
  FALSE
}

read_onecol_csv <- function(path) {
  v <- suppressWarnings(vroom::vroom(path, col_names = FALSE, delim = ","))
  as.character(v[[1]])
}

finite_range <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NULL)
  c(min(x), max(x))
}

# NEW: helper to guess column names by common candidates (case-insensitive)
guess_col <- function(cols, candidates) {
  if (!length(cols)) return(NULL)
  cols_l <- tolower(cols)

  # exact match first (in order of candidates)
  for (cand in candidates) {
    j <- which(cols_l == tolower(cand))
    if (length(j)) return(cols[j[1]])
  }
  # then "contains" match
  for (cand in candidates) {
    j <- which(grepl(tolower(cand), cols_l, fixed = TRUE))
    if (length(j)) return(cols[j[1]])
  }
  cols[1]
}

safe_ttest_p <- function(x, g, paired = FALSE, var.equal = FALSE) {
  x <- as.numeric(x)
  g <- as.factor(g)
  if (length(unique(g)) != 2) return(1)
  a <- x[g == levels(g)[1]]
  b <- x[g == levels(g)[2]]
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) return(1)
  if (stats::sd(a) == 0 && stats::sd(b) == 0) return(1)
  out <- try(stats::t.test(x ~ g, paired = paired, var.equal = var.equal)$p.value, silent = TRUE)
  if (inherits(out, "try-error") || !is.finite(out)) 1 else as.numeric(out)
}

safe_wilcox_p <- function(x, g, paired = FALSE) {
  x <- as.numeric(x)
  g <- as.factor(g)
  if (length(unique(g)) != 2) return(1)
  a <- x[g == levels(g)[1]]
  b <- x[g == levels(g)[2]]
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 1 || length(b) < 1) return(1)
  out <- try(stats::wilcox.test(x ~ g, paired = paired)$p.value, silent = TRUE)
  if (inherits(out, "try-error") || !is.finite(out)) 1 else as.numeric(out)
}

# Parse mzMine-like feature table:
# - features are rows
# - samples are columns matching keywords
# - returns sample x feature matrix + feature map (id, mz, rt, Feature)
parse_feature_table_to_matrix <- function(raw_df,
                                         row_id_col,
                                         mz_col,
                                         rt_col,
                                         sample_keywords) {
  raw_df <- clean_mzmine_export(raw_df)
  cols <- names(raw_df)

  validate(
    need(row_id_col %in% cols, "Row ID column not found."),
    need(mz_col %in% cols,     "m/z column not found."),
    need(rt_col %in% cols,     "rt column not found.")
  )

  sidx <- multi_sample_idx(cols, sample_keywords)
  validate(need(length(sidx) > 0,
                sprintf("No sample columns matched keywords: %s",
                        paste(sample_keywords, collapse = ", "))))

  sample_cols <- cols[sidx]

  # samples x features
  mat <- as.data.frame(data.table::transpose(raw_df[, sample_cols, drop = FALSE]),
                       check.names = FALSE, stringsAsFactors = FALSE)
  rownames(mat) <- clean_sample_names(sample_cols)

  # feature definition
  id <- as.character(raw_df[[row_id_col]])
  mz <- suppressWarnings(as.numeric(raw_df[[mz_col]]))
  rt <- suppressWarnings(as.numeric(raw_df[[rt_col]]))

  feat_raw <- paste0(mz, "@", rt)
  Feature <- make.unique(feat_raw, sep = "_")
  colnames(mat) <- Feature

  fmap <- tibble(
    id = id,
    mz = mz,
    RT = rt,
    Feature = Feature
  )

  mat[] <- lapply(mat, function(z) suppressWarnings(as.numeric(z)))
  mat[is.na(mat)] <- 0

  list(mat = mat, fmap = fmap, raw = raw_df)
}

# Imputation (MVI)
impute_lod_random <- function(X,
                              noise_mode = c("quantile", "manual"),
                              noise_quantile = 0.25,
                              noise_manual = 50,
                              sd_val = 30,
                              seed = 1234) {
  noise_mode <- match.arg(noise_mode)
  X <- as.matrix(X)
  X[X == 0] <- NA

  nz <- 1:min(X, na.rm = T)
  if (!length(nz)) return(X * 0)

  if (noise_mode == "quantile") {
    noise <- as.numeric(stats::quantile(nz, probs = noise_quantile, na.rm = TRUE))
  } else {
    noise <- as.numeric(noise_manual)
  }
  noise <- ifelse(is.finite(noise) && noise > 0, noise, 1)

  sd_val <- as.numeric(sd_val)
  sd_val <- ifelse(is.finite(sd_val) && sd_val >= 0, sd_val, 0)

  set.seed(seed)
  X[is.na(X)] <- 0
  idx <- X == 0
  imp <- abs(stats::rnorm(sum(idx), mean = noise, sd = sd_val))
  X[idx] <- imp
  X
}

# Stats + FC into long volcano table
compute_stats_long <- function(df_used,
                               test = c("Student", "Wilcoxon"),
                               adj  = c("BH","holm","hochberg","hommel","bonferroni","BY","fdr","none"),
                               paired = FALSE,
                               eqvar  = FALSE,
                               pseudocount = 1.1) {
  test <- match.arg(test)
  adj  <- match.arg(adj)

  validate(need("Label" %in% names(df_used), "Internal error: Label column missing."))

  feats <- setdiff(colnames(df_used), "Label")
  validate(need(length(feats) > 0, "No feature columns detected."))

  gr <- as.factor(df_used$Label)
  lev <- levels(gr)
  validate(need(length(lev) >= 2, "Need at least 2 Label groups to run statistics."))

  comb <- t(combn(lev, 2)) # rows: (A,B) in level order
  out_list <- vector("list", nrow(comb))

  for (i in seq_len(nrow(comb))) {
    gden <- comb[i, 1]  # denominator
    gnum <- comb[i, 2]  # numerator
    comp <- paste0(gnum, " / ", gden)

    sub <- df_used[df_used$Label %in% c(gden, gnum), c("Label", feats), drop = FALSE]
    sub$Label <- factor(as.character(sub$Label), levels = c(gden, gnum))

    # p-values on RAW scale (never log)
    if (test == "Student") {
      p <- vapply(feats, function(f) safe_ttest_p(sub[[f]], sub$Label, paired = paired, var.equal = eqvar), numeric(1))
    } else {
      p <- vapply(feats, function(f) safe_wilcox_p(sub[[f]], sub$Label, paired = paired), numeric(1))
    }
    padj <- if (adj == "none") p else stats::p.adjust(p, method = adj)

    # group means on RAW scale
    Xnum <- sub[sub$Label == gnum, feats, drop = FALSE]
    Xden <- sub[sub$Label == gden, feats, drop = FALSE]
    mean_num_raw <- colMeans(Xnum, na.rm = TRUE)
    mean_den_raw <- colMeans(Xden, na.rm = TRUE)

    # FC 
    mean_num_log2 <- log2(colMeans(as.matrix(Xnum), na.rm = TRUE) + pseudocount)
    mean_den_log2 <- log2(colMeans(as.matrix(Xden), na.rm = TRUE) + pseudocount)
    FC <- mean_num_log2 - mean_den_log2

    Mean <- 0.5 * (mean_num_raw + mean_den_raw)

    dd <- tibble(
      Groups = comp,
      Group_num = gnum,
      Group_den = gden,
      Feature = feats,
      `Adj.p-value` = as.numeric(padj),
      Mean = as.numeric(Mean),
      mean_num = as.numeric(mean_num_raw),
      mean_den = as.numeric(mean_den_raw),
      FC = as.numeric(FC)
    )
    dd$`Adj.p-value.log` <- -log10(pmax(dd$`Adj.p-value`, .Machine$double.xmin))
    dd$Significant <- (dd$`Adj.p-value` <= 0.05) & (abs(dd$FC) >= 1)

    out_list[[i]] <- dd
  }

  dplyr::bind_rows(out_list)
}

# ----------------------------- UI -----------------------------------------

ui <- fluidPage(
  useShinyjs(),

tags$head(tags$style(HTML("
  .app-footer { position: fixed; left:0; right:0; bottom:0; 
                text-align:center; font-size:12px; opacity:0.75;
                padding:8px; background: rgba(255,255,255,0.8);
                border-top: 1px solid #ddd; z-index: 9999; }
  body { padding-bottom: 45px; }
"))),

tags$head(
  tags$title("Metabocano — Volcano Explorer"),
  tags$link(rel = "icon", type = "image/png",
            href = "https://raw.githubusercontent.com/plyush1993/Metabocano/main/sticker.png")
),

div(
  class = "app-footer",
  HTML('Created by: Ivan Plyushchenko &nbsp;|&nbsp;
       <a href="https://github.com/plyush1993/Metabocano" target="_blank">GitHub repository</a>')
),
  
  div(
  style = "
    width: 100%;
    display: flex;
    align-items: center;
    justify-content: center;
    margin-bottom: 20px;
  ",
  
  tags$img(
    src = 'https://raw.githubusercontent.com/plyush1993/Metabocano/main/sticker.png',
    height = '150px',
    style = 'margin-right: 20px;'
  ),
  
  div(
    style = '
      font-size: 32px;
      font-weight: 900;
      color: #66CDAA; 
      text-align: center;
    ',
    "Enhanced Interactive Volcano Plot for Metabolomics Studies"
  )
), #c("#48D1CC", "#66CDAA", "#00CDCD")
  
 theme = shinytheme("flatly"), 
  setBackgroundColor(color = c("#FFFFFF", "#FFFFFF", "#67CFAC61"), gradient = "linear", direction = "bottom"),

  tags$head(
    tags$style(HTML("
      .nav-tabs > li > a {
        font-size: 20px !important;
        font-weight: bold !important;
        padding: 12px 18px !important;
      }
      .nav-tabs > li.active > a {
        font-size: 22px !important;
      }
    "))
  ),

  tags$head(tags$style(HTML("
    .shiny-output-error-validation { color:#b00020 !important; font-size: 18px !important; font-weight:800 !important; padding:10px; }
    .highlight { background:#fff; border:2px solid #000; color:#000; padding:8px; font-size: 18px; border-radius:8px; font-weight:bold; }
    .nav-tabs>li>a { font-size: 18px; padding: 10px 14px; }
    .small-note { font-size: 13px; opacity: 0.85; }
  "))),

  tabsetPanel(id = "tabs",
    tabPanel("1) Load & Process", value = "load",
      sidebarLayout(
        sidebarPanel(
          h3(class = "highlight", "Upload"),
          fileInput("file_data", "Upload feature table (.csv)", accept = ".csv"),

          # NEW: software tool selector (default = mzMine behavior)
          selectInput(
            "software_tool",
            "Software tool:",
            choices = c("mzMine" = "mzmine", "xcms" = "xcms", "MS-DIAL" = "msdial"),
            selected = "mzmine"
          ),

          tags$hr(),
          uiOutput("col_pickers"),

          h4("Sample columns detection"),
          selectizeInput(
            "sample_keywords",
            "Sample column keywords (pick/add multiple):",
            choices  = c(".mzML", ".mzXML", ".raw", " Peak area"),
            selected = c(".mzML", ".mzXML"),
            multiple = TRUE,
            options  = list(create = TRUE, createOnBlur = TRUE,
                            placeholder = "Type to add keyword and press Enter")
          ),

          tags$hr(),
          h3(class = "highlight", "Labels"),
          radioButtons(
            "label_source",
            "Label source:",
            choices = c("From sample names (token)" = "token",
                        "From uploaded labels CSV (1 column, no header)" = "csv"),
            selected = "token"
          ),
          conditionalPanel(
            condition = "input.label_source == 'token'",
            textInput("token_sep", "Token separator", value = "_"),
            numericInput("token_index", "Token index (1-based)", value = 2, min = 1, step = 1),
            checkboxInput("clean_sample_names", "Clean sample names (remove extension/Peak area)", TRUE)
          ),
          conditionalPanel(
            condition = "input.label_source == 'csv'",
            fileInput("file_labels", "Upload labels CSV", accept = ".csv")
          ),
          checkboxInput("show_labels_table", "Show labels table", TRUE),
          tags$hr(),

          h3(class = "highlight", "Annotation (optional)"),
          materialSwitch("use_sirius", "Join SIRIUS CANOPUS summary", value = FALSE, status = "success"),
          conditionalPanel(
            condition = "input.use_sirius",
            fileInput("file_sirius", "Upload canopus_formula_summary.csv", accept = ".csv"),
            uiOutput("sirius_pickers")
          ),

          tags$hr(),
          h3(class = "highlight", "MVI (imputation)"),
          radioButtons("do_mvi", "Imputation:", c("No"="no", "Yes"="yes"), selected = "yes", inline = TRUE),
          conditionalPanel(
            condition = "input.do_mvi == 'yes'",
            radioButtons("mvi_type", "Type:", c("LOD random"="lod", "BPCA"="bpca"), selected = "lod"),
            conditionalPanel(
              condition = "input.mvi_type == 'lod'",
              radioButtons("noise_mode", "Noise:",
                           c("Quantile of non-zero values"="quantile", "Manual value"="manual"),
                           selected = "quantile"),
              conditionalPanel(condition = "input.noise_mode == 'quantile'",
                               numericInput("noise_quantile", "Quantile (0-1)", value = 0.25, min = 0, max = 1, step = 0.01)),
              conditionalPanel(condition = "input.noise_mode == 'manual'",
                               numericInput("noise_manual", "Noise value", value = 50, min = 0, step = 1)),
              numericInput("noise_sd", "SD for random values", value = 30, min = 0, step = 1)
            )
          ),

          tags$hr(),
          h3(class = "highlight", "Statistics"),
          selectInput("test_type", "Test:", c("Student", "Wilcoxon"), selected = "Student"),
          selectInput("p_adjust", "p-adjust:", c("BH","holm","hochberg","hommel","bonferroni","BY","fdr","none"), selected = "BH"),
          checkboxInput("paired", "Paired test", FALSE),
          conditionalPanel(
            condition = "input.test_type == 'Student'",
            checkboxInput("eqvar", "Equal variances (Student t-test)", FALSE)
          ),

          tags$hr(),
          actionButton("run_proc", "Run preprocessing", class = "btn btn-success"),
          tags$br(), tags$br(),
          downloadButton("dl_volcano", "Download volcano table", class = "btn-info"),
          tags$br(),tags$br(),
          downloadButton("dl_matrix", "Download processed table", class = "btn-info")
        ),

        mainPanel(
          uiOutput("raw_header"),
          DTOutput("raw_preview"),
          tags$hr(),
          uiOutput("labels_header"),
          conditionalPanel(condition = "input.show_labels_table", DTOutput("labels_table")),
          tags$hr(),
          uiOutput("proc_summary")
        )
      )
    ),

    tabPanel("2) Volcano explorer", value = "volcano",
      sidebarLayout(
        sidebarPanel(uiOutput("volcano_sidebar")),
        mainPanel(uiOutput("volcano_main"))
      )
    )
  )
)

# ----------------------------- Server -------------------------------------

server <- function(input, output, session) {

  rv <- reactiveValues(
    raw = NULL,
    mat = NULL,
    fmap = NULL,
    labels = NULL,
    df_used = NULL,   # matrix used for stats + plotting (raw or imputed)
    volcano = NULL
  )

  procReady <- reactive({
    !is.null(rv$volcano) && nrow(rv$volcano) > 0
  })

  # ---- Load raw data
raw_df <- reactive({
  req(input$file_data)
  ext <- tools::file_ext(input$file_data$name)
  validate(need(tolower(ext) == "csv", "Please upload a .csv file"))
  clean_mzmine_export(vroom::vroom(input$file_data$datapath, delim = ","))
}) %>% bindCache(input$file_data$name)

  # column pickers for raw (NOW depends on software_tool for default selections)
  output$col_pickers <- renderUI({
    req(raw_df())
    cols <- names(raw_df())

    tool <- input$software_tool %||% "mzmine"

    # candidates per tool (defaults only; user can override)
    rid_cand <- switch(tool,
      xcms   = c("feature", "featureid", "feature_id", "id", "row id", "rowID"),
      msdial = c("alignment id", "alignmentid", "peak id", "peakid", "id", "row id"),
      c("row id", "rowID", "id")
    )
    mz_cand <- switch(tool,
      xcms   = c("mzmed", "mz", "m/z", "mzmin", "mzmax"),
      msdial = c("average mz", "average m/z", "mz", "m/z"),
      c("row m/z", "row mz", "m/z", "mz")
    )
    rt_cand <- switch(tool,
      xcms   = c("rtmed", "rt", "retention time", "rtmin", "rtmax"),
      msdial = c("average rt(min)", "average rt (min)", "rt (min)", "rt", "retention time"),
      c("row retention time", "retention time", "rt")
    )

    default_rid <- guess_col(cols, rid_cand)
    default_mz  <- guess_col(cols, mz_cand)
    default_rt  <- guess_col(cols, rt_cand)

    tagList(
      selectInput("row_id_col", "Row ID column:", choices = cols, selected = default_rid %||% cols[1]),
      selectInput("mz_col",     "m/z column:",   choices = cols, selected = default_mz  %||% cols[1]),
      selectInput("rt_col",     "rt column:",    choices = cols, selected = default_rt  %||% cols[1])
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

  # ---- Build matrix + fmap
  built <- reactive({
    req(raw_df(), input$row_id_col, input$mz_col, input$rt_col, input$sample_keywords)
    parse_feature_table_to_matrix(
      raw_df(),
      row_id_col = input$row_id_col,
      mz_col     = input$mz_col,
      rt_col     = input$rt_col,
      sample_keywords = input$sample_keywords
    )
  })

  sample_names <- reactive({
    req(built())
    rownames(built()$mat)
  })

  # ---- Labels
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

  # ---- SIRIUS pickers
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
      selectInput("sirius_idcol", "SIRIUS mappingFeatureId column:",
                  choices = cols, selected = if ("mappingFeatureId" %in% cols) "mappingFeatureId" else cols[1]),
      selectInput("sirius_npcol", "NPC class column:",
                  choices = cols, selected = if ("NPC#class" %in% cols) "NPC#class" else cols[1]),
      selectInput("sirius_cfcol", "ClassyFire class column:",
                  choices = cols, selected = if ("ClassyFire#class" %in% cols) "ClassyFire#class" else cols[1])
    )
  })

  # ---- Process button
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
        if (identical(input$mvi_type, "lod")) {
          Xm <- impute_lod_random(
            X0,
            noise_mode     = input$noise_mode %||% "quantile",
            noise_quantile = input$noise_quantile %||% 0.25,
            noise_manual   = input$noise_manual %||% 50,
            sd_val         = input$noise_sd %||% 30,
            seed           = 1234
          )
        } else {
          X0m <- as.matrix(X0)
          X0m[X0m == 0] <- NA
          Xm <- MsCoreUtils::impute_matrix(t(X0m), method = "bpca")
          Xm <- t(Xm)
        }
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
        pseudocount = 1.1
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
      #tags$hr(),  
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
      plotlyOutput("volcano_plot", height = "520px"),
      div(style = "height:8px;"),
      plotlyOutput("feature_plot", height = "260px")
    )
  })

  output$volcano_sliders <- renderUI({
    req(procReady())
    dd <- rv$volcano

    mzr <- finite_range(dd$mz)
    rtr <- finite_range(dd$RT)
    mr  <- finite_range(log10(dd$Mean))
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
  
    # (optional) significant + feature selection first
    if (isTRUE(input$sig_only)) {
      dd <- dd %>% dplyr::filter(Significant)
    }
    if (!is.null(input$sel_feat) && length(input$sel_feat) > 0) {
      dd <- dd %>% dplyr::filter(Feature %in% input$sel_feat)
    }
  
    # IMPORTANT: ASSIGN the main filters back to dd
    dd <- dd %>%
      dplyr::filter(
        is.finite(mz), is.finite(RT), is.finite(Mean),
        mz >= input$mz_range[1], mz <= input$mz_range[2],
        RT >= input$rt_range[1], RT <= input$rt_range[2],
        log10(Mean + 1.1) >= input$intensity_range[1],
        log10(Mean + 1.1) <= input$intensity_range[2]
      )
  
    # FDR filter
    if (isTRUE(input$use_fdr_filter)) {
      thr <- input$fdr_cutoff %||% 1.30103  # default = -log10(0.05)
      dd <- dd %>% dplyr::filter(is.finite(`Adj.p-value.log`), `Adj.p-value.log` >= thr)
    }
  
    output$fdr_equiv_text <- renderUI({
      req(input$fdr_cutoff)
      tags$small(sprintf("Equivalent Adj.p-value ≤ %.3g", 10^(-input$fdr_cutoff)))
    })
    
    # FC filter
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
      "<br>Feature: ", dd$Feature,
      "<br>ID: ", dd$id,
      "<br>m/z: ", dd$mz,
      "<br>RT: ", dd$RT,
      "<br>NPC: ", dd$`NPC#class`,
      "<br>ClassyFire: ", dd$`ClassyFire#class`
    )

    fc_line <- if (isTRUE(input$use_fc_filter)) input$fc_thr else 1
    ythr <- if (isTRUE(input$use_fdr_filter)) (input$fdr_cutoff %||% 1.30103) else 1.30103
    
    shapes <- list()
    
    # FC lines
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
      # default volcano lines at +/-1
      shapes <- c(shapes, list(
        list(type="line", x0=-1, x1=-1, xref="x", y0=0, y1=1, yref="paper", line=list(dash="dot")),
        list(type="line", x0= 1, x1= 1, xref="x", y0=0, y1=1, yref="paper", line=list(dash="dot"))
      ))
    }
    
    # FDR/p line
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
    req(procReady(), rv$df_used)

    click <- event_data("plotly_click", source = "volcano")
    if (is.null(click) || is.null(click$key)) return(NULL)

    key <- click$key[[1]]
    parts <- str_split_fixed(key, "__", 2)
    comp  <- parts[1, 1]
    feat  <- parts[1, 2]

    df_used <- rv$df_used
    validate(need(feat %in% colnames(df_used), "Clicked feature not found in matrix."))

    row <- rv$volcano %>% filter(Groups == comp, Feature == feat) %>% slice(1)

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
        colors = "Dark2",
        type = "scatter", mode = "markers",
        color = xx,
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

shinyApp(ui, server)
#.....................................................