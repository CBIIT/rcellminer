[![Travis-CI Build Status](https://travis-ci.org/cannin/rcellminer.svg?branch=master)](https://travis-ci.org/cannin/rcellminer)
[![codecov.io](https://codecov.io/github/cannin/rcellminer/coverage.svg?branch=master)](https://codecov.io/github/cannin/rcellminer?branch=master)

# RCellminer

This R package provide access to the data and functions to analyze data from [CellMiner](http://discover.nci.nih.gov/cellminer).

# Installation

## Install from GitHub
    setRepositories(ind=1:6)
    options(repos="http://cran.rstudio.com/")

    if (!require("devtools")) install.packages("devtools")

    source("http://bioconductor.org/biocLite.R")

    if (!require("Biobase")) {
        biocLite()
    }

    biocLite("BiocStyle")

    library(devtools)

    # Optional for interactive plots
    install_github('rCharts', 'ramnathv')

    install_github("cannin/rcellminer",
                      build_vignette=FALSE,
                      dependencies=TRUE,
                      args="--no-multiarch",
                      ref="master")

    install_github("github/rcellminerData",
                      build_vignette=FALSE,
                      dependencies=TRUE,
                      args="--no-multiarch",
                      ref="master")

## Install from local
### Open Project

Open .Rproj file in RStudio

### Install Dependencies
    if (require("devtools")) install.packages("devtools")

    setRepositories(ind=1:6)
    library(devtools)

    install_deps(".")

# Run Vignette

Open vignette (Rmd files) from vignettes/ folder and use "Knit HTML" button to generate the HTML file.
