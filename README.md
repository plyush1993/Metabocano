[![Project Status:](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://img.shields.io/badge/R≥4.5.0-5fb9ed.svg?style=flat&logo=r&logoColor=white?)](https://cran.r-project.org/index.html)
[![License](https://img.shields.io/badge/GPLv3-indianred.svg?style=flat&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
# Metabocano <img src="inst/www/sticker.png" align="right" height="180" width="160">

### Description :bookmark_tabs:
The [`Shiny App`](https://shiny.posit.co/) for making an enhanced interactive volcano plot for metabolomics studies.
- Directly reads the output peak table from [`mzMine`](https://mzio.io/mzmine-news/), [`xcms`](https://www.bioconductor.org/packages/release/bioc/html/xcms.html), [`MS-DIAL`](https://systemsomicslab.github.io/compms/msdial/main.html), and Default format (see [`examples of inputs`](https://github.com/plyush1993/Metabocano/tree/main/toy_examples))
- Merges with the annotation results from [`SIRIUS`](https://bio.informatik.uni-jena.de/software/sirius/), [`GNPS`](https://gnps2.org/homepage)
- Intersects peak table with the annotation table from [`SIRIUS`](https://bio.informatik.uni-jena.de/software/sirius/), and Pairs List from [`GNPS`](https://gnps2.org/homepage)
- Performs imputation by noise, and statistical tests by group pairs
- Generates several outputs:
  - interactive [`Plotly`](https://plotly.com/)-type volcano plot with built-in filters, sliders, and formatting
  - table reformatted for [`MetaboAnalyst`](https://www.metaboanalyst.ca/home.xhtml) input
  - zip archive for [`AutoPlotter`](https://mpietzke.shinyapps.io/AutoPlotter/) *Compounds in Columns* input
  - table with statistics and [`SIRIUS`](https://bio.informatik.uni-jena.de/software/sirius/)/[`GNPS`](https://gnps2.org/homepage) annotations, which can be merged in [`Cytoscape`](https://cytoscape.org/) by *id* column comes from the used peak table

### Launch the App :rocket:
**Shiny deployment**<br>
[**`https://plyush1993.shinyapps.io/Metabocano/`**](https://plyush1993.shinyapps.io/Metabocano) <br><br>
**Run locally**<br>
Install:
```r
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("plyush1993/metabocano", INSTALL_opts = "--no-multiarch")
```
or
```r
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pak("plyush1993/metabocano")
```
Run:
```r
metabocano::run_metabocano()
```
<br>

> [!IMPORTANT]
>The App was compiled using [`R version 4.5.0`](https://cran.r-project.org/bin/windows/base/old/4.5.0/) 
<br>

### Contact :mailbox_with_mail:
Please send any comment, suggestion or question you may have to the author (Dr. Ivan Plyushchenko):  
<div> 
  <a href="mailto:plyushchenko.ivan@gmail.com"><img src="https://img.shields.io/badge/-4a9edc?style=for-the-badge&logo=gmail" height="28" alt="Email" /></a>
  <a href="https://github.com/plyush1993"><img src="https://img.shields.io/static/v1?style=for-the-badge&message=%20&color=181717&logo=GitHub&logoColor=FFFFFF&label=" height="28" alt="GH" /></a>
  <a href="https://orcid.org/0000-0003-3883-4695"><img src="https://img.shields.io/badge/-A6CE39?style=for-the-badge&logo=ORCID&logoColor=white" height="28" alt="ORCID" /></a>
</div>

