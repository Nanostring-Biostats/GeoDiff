# GeoDiff

## Overview

We present GeoDiff, an R package for count generating models for analyzing Geomx RNA data, demonstrated by example data analysis pipeline vignette. The package provides a series of statistical models using count generating distributions for fundamental tasks including background modelling, target and sample QC, normalization and differential expression analysis. The new methods are shown to be more suitable for these tasks in real examples than some conventional methods. 


It contains the definition of the NanoStringGeoMxSet class which
inherits from Biobase’s ExpressionSet class and NanoStringRCCSet class.

## Installation

### Install the development version from GitHub
``` r
install.packages("devtools")
library("devtools")
devtools::install_github("Nanostring-Biostats/GeoDiff", ref = "dev")
```

## License
This project is licensed under the [MIT license](LICENSE).
