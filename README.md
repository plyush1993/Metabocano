# Metabocano <img src="sticker.png" align="right" height="180" width="160">

### Description :bookmark_tabs:
The [Shiny App](https://shiny.posit.co/) for making an enhanced interactive volcano plot for metabolomics studies.
The App is available to directly read the output peak table from [mzMine](https://mzio.io/mzmine-news/), [xcms](https://www.bioconductor.org/packages/release/bioc/html/xcms.html), [MS-DIAL](https://systemsomicslab.github.io/compms/msdial/main.html), and Default (see [examples of inputs](https://github.com/plyush1993/Metabocano/tree/main/toy_examples)) then join it with the annotation table from [SIRIUS](https://bio.informatik.uni-jena.de/software/sirius/), perform imputation, and statistical tests. After preprocessing, the App generates several outputs: interactive [Plotly](https://plotly.com/)-type volcano plot, the table reformatted for [MetaboAnalyst](https://www.metaboanalyst.ca/home.xhtml) input, and a statistical table with SIRIUS annotation.

### Launch the App :rocket:
Shiny deployment:<br>
[**`https://plyush1993.shinyapps.io/Metabocano/`**](https://plyush1993.shinyapps.io/Metabocano) <br><br>
Run locally:
```r
cat("Checking required packages (auto-installing if missing)\n")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load("shiny", "DT", "shinythemes", "shinyWidgets", "shinyjs", "vroom", "dplyr", "data.table", "tidyr", "stringr", "plotly", "viridisLite", "tibble", "shinyBS")

source("https://raw.githubusercontent.com/plyush1993/Metabocano/refs/heads/main/app.R")
shiny::shinyApp(ui, server)
```
<br>

> [!IMPORTANT]
>The [App's script](https://github.com/plyush1993/Metabocano/blob/main/app.R) was compiled using [R version 4.1.2](https://cran.r-project.org/bin/windows/base/old/4.1.2/) 
<br>

### Contact :mailbox_with_mail:
Please send any comment, suggestion or question you may have to the author (Dr. Ivan Plyushchenko):  
<div> 
  <a href="mailto:plyushchenko.ivan@gmail.com"><img src="https://img.shields.io/badge/-4a9edc?style=for-the-badge&logo=gmail" height="28" alt="Email" /></a>
  <a href="https://github.com/plyush1993"><img src="https://img.shields.io/static/v1?style=for-the-badge&message=%20&color=181717&logo=GitHub&logoColor=FFFFFF&label=" height="28" alt="GH" /></a>
  <a href="https://orcid.org/0000-0003-3883-4695"><img src="https://img.shields.io/badge/-A6CE39?style=for-the-badge&logo=ORCID&logoColor=white" height="28" alt="ORCID" /></a>
</div>

