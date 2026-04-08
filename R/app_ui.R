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

tags$head(
  tags$title("Metabocano"),
  tags$link(rel = "icon", type = "image/png",
            href = "www/sticker.png")
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
    src = 'www/sticker.png',
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
          tags$hr(),
          uiOutput("col_pickers"),

          h4("Sample columns detection"),
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
          class = "",
          style = "font-weight: bold; margin-left: 10px; margin-top: -2px;"
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
          selectInput("test_type", "Test:", c("Student", "Wilcoxon"), selected = "Student"),
          selectInput("p_adjust", "p-adjust:", c("BH","holm","hochberg","hommel","bonferroni","BY","fdr","none"), selected = "BH"),
          checkboxInput("paired", "Paired test", FALSE),
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
          tags$hr(),
          actionButton("run_proc", "Run preprocessing", class = "btn btn-success"),
          tags$br(), tags$br(),
          downloadButton("dl_volcano", "Download volcano table", class = "btn-info"),
          actionButton("btn1", "?"),
          bsTooltip("btn1",
          title = "<b>Download table with all calculated statistical values.</b>", "right", trigger = "click", options = list(container = "body")),

          tags$br(),tags$br(),
          downloadButton("dl_matrix", "Download processed table", class = "btn-info"),
          actionButton("btn2", "?"),
          bsTooltip("btn2",
          title = "<b>Download peak table after all processing steps and <em>Label</em> column.</b><br>Suitable as input in MetaboAnalyst (www.metaboanalyst.ca/).", "right", trigger = "click", options = list(container = "body"))
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
}
