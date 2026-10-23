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
#' @import igraph
#' @import ComplexHeatmap
#' @import InteractiveComplexHeatmap
#' @import ggrepel
app_ui <- function() {
fluidPage(
  useShinyjs(),

  tags$head(
  tags$script(HTML("
    (function () {

      function fixFileInputs() {
        document.querySelectorAll('input[type=file]').forEach(
          function (input) {

            const button = input.closest('.btn-file');

            if (!button) return;

            // Keep the hidden input inside its Browse button.
            button.style.setProperty(
              'position', 'relative', 'important'
            );

            input.style.setProperty(
              'position', 'absolute', 'important'
            );
            input.style.setProperty(
              'top', '0', 'important'
            );
            input.style.setProperty(
              'left', '0', 'important'
            );
            input.style.setProperty(
              'width', '1px', 'important'
            );
            input.style.setProperty(
              'height', '1px', 'important'
            );
            input.style.setProperty(
              'opacity', '0', 'important'
            );
            input.style.setProperty(
              'pointer-events', 'none', 'important'
            );
          }
        );
      }

      // Inputs present when the page first loads.
      $(fixFileInputs);

      // Inputs subsequently created by renderUI().
      $(document).on('shiny:bound', fixFileInputs);

    })();
  "))
),

tags$head(tags$style(HTML("
  .app-footer { position: fixed; left:0; right:0; bottom:0;
                text-align:center; font-size:12px; opacity:0.75;
                padding:8px; background: rgba(255,255,255,0.8);
                border-top: 1px solid #ddd; z-index: 9999; }
  body { padding-bottom: 45px; }
"))),

tags$head(
  tags$script(
    HTML("
      $(document).on(
        'click',
        '#selected_feature_info_btn',
        function() {

          var txt = $(this).attr('data-copy');

          if (!txt || !navigator.clipboard) {
            return;
          }

          navigator.clipboard.writeText(txt).then(function() {
            Shiny.setInputValue(
              'selected_feature_info_copied',
              Date.now(),
              {priority: 'event'}
            );
          });
        }
      );
    ")
  )
),

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

div(class = "app-footer", HTML('
    <span class="footer-text">by Plyushchenko I.V.</span>
    <span class="footer-sep">&nbsp;|&nbsp;</span>
    <span class="footer-text">GPLv3</span>
    <span class="footer-sep">&nbsp;|&nbsp;</span>
     <a id="latest-release-link"
     class="footer-link"
     href="https://github.com/plyush1993/metabocano/releases/latest"
     target="_blank">v. </a>

    <script>
    fetch("https://api.github.com/repos/plyush1993/metabocano/releases/latest")
      .then(function(response) {
        if (!response.ok) throw new Error("GitHub release request failed");
        return response.json();
      })
      .then(function(data) {
        var link = document.getElementById("latest-release-link");
        if (link && data.tag_name) {
          link.textContent = "v. " + data.tag_name;
          if (data.html_url) {
            link.href = data.html_url;
          }
        }
      })
      .catch(function(error) {
        var link = document.getElementById("latest-release-link");
        if (link) {
          link.textContent = "Latest release";
          link.href = "https://github.com/plyush1993/metabocano/releases/latest";
        }
      });
  </script>

  ')),

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
  choices = c(
    "mzMine" = "mzmine",
    "xcms" = "xcms",
    "MS-DIAL" = "msdial",
    "Generic" = "default"
  ),
  selected = "mzmine"
),

fileInput(
  "file_data",
  "Upload feature table (.csv)",
  accept = ".csv"
),

helpText(
  HTML(
    "<i class='fa fa-info-circle'></i> Need data to test? Download example datasets from our <a href='https://github.com/plyush1993/Metabocano' target='_blank'>GitHub</a>."
  )
),

uiOutput("upload_tab_error"),

uiOutput("col_pickers"),

          radioButtons(
  "sample_mode",
  "How to define sample columns?",
  choices = c(
    "Auto-detect numeric sample columns" = "auto",
    "By keyword match" = "kws",
    "Pick columns manually" = "manual"
  ),
  selected = "kws"
),

conditionalPanel(
  condition = "input.sample_mode == 'kws'",

  selectizeInput(
    "sample_keywords",
    "Sample column keywords (pick/add multiple):",
    choices = c(
      ".mzML", ".mzXML", ".raw", ".d", ".wiff", ".lcd",
      "Peak area", "Area", "_Area"
    ),
    selected = c(".mzML", ".mzXML"),
    multiple = TRUE,
    options = list(
      create = TRUE,
      createOnBlur = TRUE,
      placeholder = "Type to add keyword and press Enter"
    )
  )
),

conditionalPanel(
  condition = "input.sample_mode == 'manual'",
  uiOutput("manual_sample_cols_ui")
),

          tags$hr(),
          h3(class = "highlight", "Labels"),
          helpText(HTML("<i class='fa fa-info-circle'></i> Avoid double underscores (`__`) in group labels")),

radioButtons(
  "label_source",
  "Label source:",
  choices = c(
    "From sample names (by token)" = "token",
    "From metadata CSV (match by sample name)" = "metadata",
    "From uploaded labels CSV (1 column, no header)" = "csv",
    "Manual editable table" = "manual"
  ),
  selected = "token"
),

conditionalPanel(
  condition = "input.label_source == 'token' || input.label_source == 'manual'",
  textInput("token_sep", "Token separator", value = "_"),
  numericInput("token_index", "Token index (1-based)", value = 2, min = 1, step = 1)
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
    value = TRUE
  ),

  conditionalPanel(
    condition = "input.metadata_clean_sample_names == true",

    selectizeInput(
      "metadata_remove_suffixes",
      "Remove suffixes/extensions:",
      choices = c(
        ".mzML", ".mzXML", ".raw", ".RAW", ".lcd",
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

          h3(class = "highlight", "Join with Annotation"),

annotation_panel(
  switch_id = "use_peak_extra_cols",
  label = "Additional peak-table columns",
  tooltip_id = "btnAD",

  tooltip_text = paste0(
    "Selected columns will be added to the downloaded volcano table ",
    "and displayed after clicking a volcano point."
  ),

  uiOutput("peak_extra_cols_ui")
),

annotation_panel(
  switch_id = "use_sirius",
  label = "SIRIUS / CANOPUS annotation",
  tooltip_id = "btn5",

  tooltip_text = paste0(
    "<b>Join SIRIUS annotations with the processed table.</b><br>",
    "The selected peak-table Feature ID column is matched ",
    "to the selected SIRIUS mapping ID column.<br>",
    "Default mapping ID: <em>mappingFeatureId</em><br>",
    "Default NPC column: <em>NPC#class</em><br>",
    "Default ClassyFire column: <em>ClassyFire#class</em>"
  ),

  fileInput(
    "file_sirius",
    "Upload SIRIUS output (.csv/.tsv/.txt)",
    accept = c(".csv", ".tsv", ".txt")
  ),

  uiOutput("sirius_pickers")
),

annotation_panel(
  switch_id = "use_gnps_annotation",
  label = "GNPS library annotation",
  tooltip_id = "btn_gnps_annotation",

  tooltip_text = paste0(
    "<b>Join GNPS library annotations with the processed table.</b><br>",
    "The selected peak-table Feature ID column is matched ",
    "to the selected GNPS ID column.<br>",
    "Default GNPS ID: <em>#Scan#</em><br>",
    "Default annotation: <em>Compound_Name</em>"
  ),

  fileInput(
    "file_gnps_annotation",
    "Upload GNPS library results (.tsv/.txt/.csv)",
    accept = c(".tsv", ".txt", ".csv")
  ),

  uiOutput("gnps_annotation_pickers")
),

annotation_panel(
  switch_id = "use_main_gnps_pairs",
  label = "GNPS network / ComponentIndex",
  tooltip_id = "btn_main_gnps_pairs",

  tooltip_text = paste0(
    "<b>Enable filtering by GNPS network component.</b><br>",
    "Select the peak-table matching ID column, both pairs-file ",
    "node ID columns, and the component column.<br>",
    "Default pairs columns: <em>CLUSTERID1</em>, ",
    "<em>CLUSTERID2</em>, and <em>ComponentIndex</em>.<br>",
    "Both endpoints are assigned to their component."
  ),

  fileInput(
    "file_main_gnps_pairs",
    "Upload GNPS network pairs (.tsv/.txt/.csv)",
    accept = c(".tsv", ".txt", ".csv")
  ),

  uiOutput("main_gnps_pairs_pickers")
),

annotation_panel(
  switch_id = "use_other_annotation",
  label = "Other annotation source",
  tooltip_id = "btn_other_annotation",

  tooltip_text = paste0(
    "<b>Join annotations from an external table.</b><br>",
    "Choose a peak-table ID column and the corresponding ",
    "ID column in the annotation file.<br>",
    "Choose one primary annotation column and optionally ",
    "add additional columns."
  ),

  fileInput(
    "file_other_annotation",
    "Upload annotation table (.csv/.tsv/.txt)",
    accept = c(".csv", ".tsv", ".txt")
  ),

  uiOutput("other_annotation_pickers")
),

          tags$hr(),
          h3(class = "highlight", "Imputation by Noise"),
          radioButtons("do_mvi", "Imputation:", c("No"="no", "Yes"="yes"), selected = "no", inline = TRUE),
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
          radioButtons("comparison_mode", "Comparison selection:",
            choices = c(
              "Reference group vs all others" = "reference",
              "Choose comparisons manually" = "manual"), selected = "reference"),
          uiOutput("comparison_picker"),
          selectInput(
            "test_type", "Test:",
            c("Student", "Wilcoxon", "limma (Moderated t-test)" = "limma"), selected = "Student"),
          selectInput("p_adjust", "p-adjust:", c("BH","holm","hochberg","hommel","bonferroni","BY","fdr","none"), selected = "BH"),
          conditionalPanel(
            condition = "input.test_type == 'Student' || input.test_type == 'Wilcoxon'",
            checkboxInput(
  "paired",
  "Paired test — samples paired by order",
  FALSE
),

conditionalPanel(
  condition = "input.paired == true",

  helpText(
    paste(
      "Samples are paired by their order within each group.",
      "Check every pair below before running preprocessing.",
      "Sample names are not used to identify matching subjects."
    )
  ),

  tags$details(
    tags$summary("Show sample pairs"),
    div(
      style = "overflow-x: auto;",
      tableOutput("paired_sample_preview")
    )
  )
)
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
          downloadButton("dl_annotation", "Feature table csv", class = "btn-info"),
          actionButton("btn_annotation", "?"),
          bsTooltip("btn_annotation",
            title = paste0(
              "<b>Download table with one row per feature with annotation information.</b><br>",
              "Always includes <em>Feature ID</em>, <em>Annotation matching ID</em>, ",
              "<em>m/z</em>, and <em>RT</em>.<br>",
              "After preprocessing, all joined Peak table, SIRIUS, GNPS, ",
              "and Other Annotation columns are also included."
            ),
            placement = "right",
            trigger = "click",
            options = list(container = "body")
          ),
          tags$br(), tags$br(),
          downloadButton("dl_volcano", "Volcano table csv", class = "btn-info"),
          actionButton("btn1", "?"),
          bsTooltip("btn1",
          title = "<b>Download table with all calculated statistical values. Can be merged <em>GNPS-derived .cys file</em> in Cytoscape by <em>id</em> column.</b>", "right", trigger = "click", options = list(container = "body")),

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
          uiOutput("label_upload_warning"),
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
        mainPanel(
  conditionalPanel(
    condition = "output.volcano_ready === 'yes'",
    volcano_main_ui()
  )
)
      )
    ),

navbarMenu(
  title = "3) Other Utils",

  tabPanel(
    title = "SIRIUS & GNPS merging",
    value = "sirius_gnps",

    sidebarLayout(
      sidebarPanel(
        uiOutput("sirius_gnps_sidebar")
      ),
      mainPanel(
        uiOutput("sirius_gnps_main")
      )
    )
  )
)

  )
)
}
