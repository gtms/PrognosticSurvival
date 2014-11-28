## This script installs all required packages for analysis
## Gil Tomás, IRIBHM, ULB
## <gil.tomas@ulb.ac.be>

## install Bioconductor
source ("http://bioconductor.org/biocLite.R")
biocLite ()

## install Bioconductor dependencies for MicroarrayToolbox
bioc.pkgs <- c ("annotate", "org.Hs.eg.db", "GO.db", "impute")
source ("http://bioconductor.org/biocLite.R")
lapply (bioc.pkgs, biocLite)

## install CRAN dependencies for MicroarrayToolbox
cran.pkgs <- c ("data.table", "MASS", "ggplot2", "RColorBrewer", "samr", "devtools", "ProjectTemplate")
lapply (cran.pkgs, install.packages, repos = "http://cran.r-project.org")

## install MicroarrayToolbox
library (devtools)
install_github ("gtms/MicroarrayToolbox")

## load MicroarrayToolbox and install remaining packages
library (MicroarrayToolbox)
pkgs.dfr <- data.frame (pkg = c ("pcaMethods", "qvalue", "SNAGEE",
                            "amap", "survival", "doMC", "MASS", "stargazer"),
                        repos = c (rep ("bioc", 3), rep ("cran", 5)),
                        stringsAsFactors = FALSE)
instant.pkgs (pkgs.dfr, verbose = TRUE)

## quit
sessionInfo ()
q (save = "no")

