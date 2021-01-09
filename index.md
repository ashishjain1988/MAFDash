MAFDash
------------------------------------------------------------
*Once you call the variants, it's a MAFDash to the finish line*


**An easy-to-use dashboard builder for mutation data**

### Installation from Github
```
install.packages(devtools)
library(devtools)
devtools::install_github("mtandon09/MAFDashRPackage")
```

### Description
[Mutation Annotation Format (MAF)](https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format/) is a tabular data format used for storing genetic mutation data. For example, [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) project has made MAF files from each project publicly available.

MAFDash is an R package that helps to quickly create an HTML dashboard to summarize and visualize data from MAF files. The resulting HTML file serves as a self-contained report that can be used to explore and share the results.

Currently, MAFDash produces mostly static plots powered by [maftools](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html),  [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) and [circlize](https://github.com/jokergoo/circlize), as well as interactive visualizations using [canvasXpress](https://cran.r-project.org/web/packages/canvasXpress/vignettes/getting_started.html) and [plotly](https://plotly.com/r/).  The report is generated with a parameterized [R Markdown](https://rmarkdown.rstudio.com/) script that uses [flexdashboard](https://rmarkdown.rstudio.com/flexdashboard/) to arrange all the information.

### How to use

- [Quick Start](articles/Quick Start.html)
- [Customizing](articles/Customizing.html)
- [Real-world example](articles/Advanced_Example.html)

### Example dashboards from the vignette
Here are some example dashboards created using TCGA data:

- [Default dashboard using TCGA-LAML](https://mtandon09.github.io/MAFDashRPackage/examples/LAML.mafdash.html)
- [Customized dashboard using TCGA-UVM](https://mtandon09.github.io/MAFDashRPackage/examples/TCGA-UVM.custom.mafdash.html.MAFDash.html)
- [Arbitrary plots using `iris` dataset](https://mtandon09.github.io/MAFDashRPackage/examples/toy_dash.html.MAFDash.html)

#### Other notes
- This repo was born out of a Shiny app, [MAFWiz](https://github.com/mtandon09/mafwiz).  Instead of relying on a Shiny server, this dashboard was an attempt to try some of those things using client-side javascript functionality.