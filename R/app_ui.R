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
app_ui <- function() {
fluidPage(
  useShinyjs(),

tags$head(tags$style(HTML("
  .app-footer { position: fixed; left:0; right:0; bottom:0;
                text-align:center; font-size:12px; opacity:0.75;
                padding:8px; background: rgba(255,255,255,0.8);
                border-top: 1px solid #ddd; z-index: 9999; }
  body { padding-bottom: 45px; }
"))),

tags$head(tags$style(HTML("
  /* Editable Labels table: prevent white-on-white editing issue */

  #labels_table.html-widget.datatables {
    background-color: transparent !important;
  }

  #labels_table .dataTables_wrapper,
  #labels_table table.dataTable,
  #labels_table .dataTables_scroll,
  #labels_table .dataTables_scrollHead,
  #labels_table .dataTables_scrollBody {
    background-color: #ffffff !important;
    color: #2c3e50 !important;
  }

  #labels_table table.dataTable th,
  #labels_table table.dataTable td {
    background-color: #ffffff !important;
    color: #2c3e50 !important;
  }

  /* Cell when focused / double-clicked / edited */
  #labels_table table.dataTable tbody td.focus,
  #labels_table table.dataTable tbody td:focus,
  #labels_table table.dataTable tbody tr.selected td,
  #labels_table table.dataTable tbody td.selected {
    background-color: #ffffff !important;
    color: #000000 !important;
    box-shadow: inset 0 0 0 2px #66CDAA !important;
  }

  /* Input box created during editing */
  #labels_table input,
  #labels_table textarea,
  #labels_table .dataTables_wrapper input,
  #labels_table .dataTables_wrapper textarea {
    background-color: #ffffff !important;
    color: #000000 !important;
    -webkit-text-fill-color: #000000 !important;
    caret-color: #000000 !important;
    border: 1px solid #66CDAA !important;
  }

  /* Keep search / info / pagination readable if shown */
  #labels_table .dataTables_length,
  #labels_table .dataTables_filter,
  #labels_table .dataTables_info,
  #labels_table .dataTables_paginate {
    color: #000000 !important;
    font-weight: bold;
    padding: 5px;
  }

  /* Colored pickerInput buttons for Volcano filters */

.bootstrap-select > .dropdown-toggle[data-id='sel_feat'] {
  background-color: #66CDAA !important;
  border-color: #45b894 !important;
  color: white !important;
  font-weight: bold !important;
}

.bootstrap-select > .dropdown-toggle[data-id='npc_filter_values'] {
  background-color: #18bc9c !important;
  border-color: #13a085 !important;
  color: white !important;
  font-weight: bold !important;
}

.bootstrap-select > .dropdown-toggle[data-id='classyfire_filter_values'] {
  background-color: #18bc9c !important;
  border-color: #13a085 !important;
  color: white !important;
  font-weight: bold !important;
}

/* Highlight selected options inside dropdown */
.bootstrap-select .dropdown-menu li.selected a {
  background-color: #dff7ef !important;
  color: #000000 !important;
  font-weight: bold !important;
}

.bootstrap-select .dropdown-menu li.selected a span.check-mark {
  color: #18bc9c !important;
}

"))),

tags$head(
  tags$title("Metabocano"),
  tags$link(rel = "icon", type = "image/png",
            href = "https://raw.githubusercontent.com/plyush1993/Metabocano/main/sticker.png")
),

tags$head(
    tags$style(HTML("
      /* Increase max-width to prevent wrapping and align text left */
      .tooltip-inner {
        max-width: none !important;
        white-space: nowrap;
        text-align: left !important;
        font-size: 18px;
      }
    "))
  ),

tags$head(tags$style(HTML("
  /* make disabled download links truly inactive */
  a.shiny-download-link.disabled,
  .shiny-download-link.disabled {
    pointer-events: none !important;
    opacity: 0.5 !important;
    cursor: not-allowed !important;
  }
"))),

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
),

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

tags$head(tags$style(HTML("
    /* Existing Footer and layout styles */
    .app-footer { position: fixed; left:0; right:0; bottom:0;
                  text-align:center; font-size:12px; opacity:0.75;
                  padding:8px; background: rgba(255,255,255,0.8);
                  border-top: 1px solid #ddd; z-index: 9999; }
    body { padding-bottom: 45px; }

    /* --- NEW: Thicker Upload Progress Bar --- */
    .progress.shiny-file-input-progress {
      height: 20px !important;
      margin-top: 10px !important;
      border-radius: 5px !important;
    }

    .progress.shiny-file-input-progress .progress-bar {
      line-height: 20px !important;
      font-size: 14px !important;
      font-weight: bold !important;
      background-color: #66CDAA !important; /* Matches your app's theme color */
    }
    /* ---------------------------------------- */

    .tooltip-inner {
      max-width: none !important;
      white-space: nowrap;
      text-align: left !important;
      font-size: 18px;
    }
  "))),

  tabsetPanel(id = "tabs",
    tabPanel("1) Load & Process", value = "load",
      sidebarLayout(
        sidebarPanel(
          h3(class = "highlight", "Upload"),
          selectInput(
            "software_tool",
            "Software tool:",
            choices = c("mzMine" = "mzmine", "xcms" = "xcms", "MS-DIAL" = "msdial", "Default" = "default"),
            selected = "mzmine"
          ),
          fileInput("file_data", "Upload feature table (.csv)", accept = ".csv"),
          helpText(HTML("<i class='fa fa-info-circle'></i> Need data to test? Download example datasets from our <a href='https://github.com/plyush1993/Metabocano' target='_blank'>GitHub</a>.")),
          uiOutput("col_pickers"),

          selectizeInput(
            "sample_keywords",
            "Sample column keywords (pick/add multiple):",
            choices  = c(".mzML", ".mzXML", ".raw", " Peak area", "Area"),
            selected = c(".mzML", ".mzXML"),
            multiple = TRUE,
            options  = list(create = TRUE, createOnBlur = TRUE,
                            placeholder = "Type to add keyword and press Enter")
          ),

          tags$hr(),
          h3(class = "highlight", "Labels"),
          helpText(HTML("<i class='fa fa-info-circle'></i> Avoid double underscores (`__`) in group labels")),

radioButtons(
  "label_source",
  "Label source:",
  choices = c(
    "From sample names (token)" = "token",
    "From metadata CSV (match by sample name)" = "metadata",
    "From uploaded labels CSV (1 column, no header)" = "csv",
    "Manual editable table" = "manual"
  ),
  selected = "token"
),

conditionalPanel(
  condition = "input.label_source == 'token' || input.label_source == 'manual'",
  textInput("token_sep", "Token separator", value = "_"),
  numericInput("token_index", "Token index (1-based)", value = 2, min = 1, step = 1),
  checkboxInput("clean_sample_names", "Clean sample names (remove extension/Peak area)", TRUE)
),

conditionalPanel(
  condition = "input.label_source == 'metadata'",

  fileInput(
    "file_metadata_labels",
    "Upload metadata CSV with column names",
    accept = ".csv"
  ),

  uiOutput("metadata_sample_col_ui"),
  uiOutput("metadata_label_col_ui"),

  checkboxInput(
    "metadata_clean_sample_names",
    "Clean sample names only for metadata matching",
    value = FALSE
  ),

  conditionalPanel(
    condition = "input.metadata_clean_sample_names == true",

    selectizeInput(
      "metadata_remove_suffixes",
      "Remove suffixes/extensions:",
      choices = c(
        ".mzML", ".mzXML", ".raw", ".RAW",
        ".cdf", ".CDF", ".mzData", ".mzdata",
        ".wiff", ".WIFF", ".d", ".D",
        " Peak area", " Peak Area",
        " Peak height", " Peak Height",
        "_Area", "_Height",
        " Area", " Height"
      ),
      selected = c(
        " Peak area", " Peak height",
        "_Area", "_Height",
        " Area", " Height"
      ),
      multiple = TRUE,
      options = list(
        create = TRUE,
        createOnBlur = TRUE,
        placeholder = "Type custom suffix and press Enter"
      )
    )
  ),

  div(
    class = "small-note",
    "Metadata rows are matched by sample name, not by row order. Cleaning is used for matching."
  )
),

conditionalPanel(
  condition = "input.label_source == 'csv'",
  fileInput("file_labels", "Upload labels CSV", accept = ".csv")
),

conditionalPanel(
  condition = "input.label_source == 'manual'",
  div(
    style = "display: inline-flex; align-items: center; gap: 6px; margin-bottom: 10px;",

    actionButton(
      "fill_manual_labels",
      label = tags$span(
        HTML("Fill editable table from<br>current token labels"),
        style = "line-height: 1.1;"
      ),
      class = "btn-success",
      style = "
        font-size: 12px;
        padding: 4px 8px;
        line-height: 1.1;
        width: 150px;
        white-space: normal;
      "
    )
  ),

  div(
    class = "small-note",
    "Double-click cells in the Label column to edit group names. The Sample column is locked."
  )
),

checkboxInput("show_labels_table", "Show labels table", TRUE),
tags$hr(),

          h3(class = "highlight", "Join with Annotation"),
      div(
        style = "display: flex; align-items: center; margin-bottom: 15px;",

        materialSwitch(
          inputId = "use_sirius",
          label = "Join SIRIUS summary",
          value = FALSE,
          status = "success",
          width = "auto"
        ),

        actionButton(
          inputId = "btn5",
          label = "?",
          class = "btn-xs",
          style = "font-weight: bold; margin-left: 10px; margin-top: -20px;"
        )
      ),
      bsTooltip(
        id = "btn5",
        title = "<b>Join SIRIUS csv output with processed table (by <em>Row ID</em> column)</b>",
        placement = "right",
        trigger = "click",
        options = list(container = "body")
      ),

      conditionalPanel(
        condition = "input.use_sirius",
        fileInput("file_sirius", "Upload SIRIUS .csv output", accept = ".csv"),
        uiOutput("sirius_pickers")
      ),

          tags$hr(),
          h3(class = "highlight", "Imputation by Noise"),
          radioButtons("do_mvi", "Imputation:", c("No"="no", "Yes"="yes"), selected = "yes", inline = TRUE),
          conditionalPanel(
            condition = "input.do_mvi == 'yes'",
            conditionalPanel(
              condition = "input.do_mvi == 'yes'",
              radioButtons("noise_mode", "Noise:",
                           c("Quantile in the range 1:min"="quantile", "Manual value"="manual"),
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
          uiOutput("ref_group_picker"),
          selectInput("test_type", "Test:", c("Student", "Wilcoxon", "limma (Moderated t-test)" = "limma"), selected = "Student"),
          selectInput("p_adjust", "p-adjust:", c("BH","holm","hochberg","hommel","bonferroni","BY","fdr","none"), selected = "BH"),
          conditionalPanel(
            condition = "input.test_type == 'Student' || input.test_type == 'Wilcoxon'",
            checkboxInput("paired", "Paired test", FALSE)
          ),
          conditionalPanel(
            condition = "input.test_type == 'Student'",
            checkboxInput("eqvar", "Equal variances (Student t-test)", FALSE)
          ),
          materialSwitch(
          "log2_test",
          "Log Transformation",
          value = FALSE,
          status = "success"
          ),
          materialSwitch(
            "standard_scaling",
            "Auto Scaling",
            value = FALSE,
            status = "success"
          ),
          tags$hr(),
          actionButton("run_proc", "Run preprocessing", class = "btn btn-success"),
          tags$br(), tags$br(),
          downloadButton("dl_volcano", "Volcano table csv", class = "btn-info"),
          actionButton("btn1", "?"),
          bsTooltip("btn1",
          title = "<b>Download table with all calculated statistical values.</b>", "right", trigger = "click", options = list(container = "body")),

          tags$br(),tags$br(),
          downloadButton("dl_matrix", "MetaboAnalyst-ready csv", class = "btn-info"),
          actionButton("btn2", "?"),
          bsTooltip("btn2",
          title = "<b>Download peak table after MVI and with <em>Label</em> column.</b><br>Suitable as input in MetaboAnalyst (www.metaboanalyst.ca/).", "right", trigger = "click", options = list(container = "body")),

        tags$br(),tags$br(),
        downloadButton(
          "dl_autoplotter_zip",
          "AutoPlotter-ready ZIP",
          class = "btn-info"
        ),
        actionButton("btn_auto", "?"),
        bsTooltip(
          "btn_auto",
          title = "<b>Download AutoPlotter-ready ZIP archive.</b><br>Suitable as input as <em>Compounds in Columns</em> in Metabolite AutoPlotter (https://mpietzke.shinyapps.io/AutoPlotter/).",
          placement = "right",
          trigger = "click",
          options = list(container = "body"))
        ),

        mainPanel(
          uiOutput("raw_header"),
          DTOutput("raw_preview"),
          tags$hr(),
          uiOutput("labels_header"),
          conditionalPanel(
  condition = "input.show_labels_table || input.label_source == 'manual'",
  DTOutput("labels_table")
),
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
}
