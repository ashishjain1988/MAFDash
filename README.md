![GitHub release (latest by date)](https://img.shields.io/github/v/release/ashishjain1988/MAFDash)
![GitHub all releases](https://img.shields.io/github/downloads/ashishjain1988/MAFDash/total)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4472978.svg)](https://doi.org/10.5281/zenodo.4472978)


MAFDash
------------------------------------------------------------
*Once you call the variants, it's a MAFDash to the finish line*


**An easy-to-use dashboard builder for mutation data**

### Installation from Github
```
# install.packages(devtools)  ## Install devtools if necessary
devtools::install_github("ashishjain1988/MAFDash")
```

### Description
[Mutation Annotation Format (MAF)](https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format/) is a tabular data format used for storing genetic mutation data. For example, [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) project has made MAF files from each project publicly available.

MAFDash is an R package that helps to quickly create an HTML dashboard to summarize and visualize data from MAF files. The resulting HTML file serves as a self-contained report that can be used to explore and share the results.

Currently, MAFDash produces mostly static plots powered by [maftools](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html),  [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) and [circlize](https://github.com/jokergoo/circlize), as well as interactive visualizations using [canvasXpress](https://cran.r-project.org/package=canvasXpress) and [plotly](https://plotly.com/r/).  The report is generated with a parameterized [R Markdown](https://rmarkdown.rstudio.com/) script that uses [flexdashboard](https://rmarkdown.rstudio.com/flexdashboard/) to arrange all the information.

### How to use

- [Quick Start](https://ashishjain1988.github.io/MAFDash/articles/Quick_Start.html)
- [Customizing](https://ashishjain1988.github.io/MAFDash//articles/Customizing.html)
- [Real-world example](https://ashishjain1988.github.io/MAFDash//articles/Advanced_Example.html)

### Example dashboards from the vignette
Here are some example dashboards created using TCGA data:

- [Default dashboard using TCGA-LAML](https://ashishjain1988.github.io/MAFDash//examples/LAML.mafdash.html)
- [Customized dashboard using TCGA-UVM](https://ashishjain1988.github.io/MAFDash//examples/TCGA-UVM.custom.mafdash.html.MAFDash.html)
- [Arbitrary plots using `iris` dataset](https://ashishjain1988.github.io/MAFDash//examples/toy_dash.html.MAFDash.html)

#### Other notes
- This repo was born out of a Shiny app, [MAFWiz](https://github.com/mtandon09/mafwiz).  Instead of relying on a Shiny server, this dashboard was an attempt to try some of those things using client-side javascript functionality.
