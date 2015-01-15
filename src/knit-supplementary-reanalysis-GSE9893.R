## This script knits the pdf for the supplementaryreanalysis of GSE9893-breast

library (knitr)

report.dir <- "reports/GSE9893-reanalysis"
setwd (report.dir)

## the root directory of the Rnw script is still the project's root directory
## opts_knit$set (root.dir = file.path (report.dir, "../.."))
opts_knit$set (root.dir = "../..")

## knit2pdf (file.path (report.dir, "supplementary-reanalysis-GSE9893.Rnw"),
##           output = file.path (report.dir, "supplementary-reanalysis-GSE9893.tex"))

knit2pdf ("supplementary-reanalysis-GSE9893.Rnw")
## output = file.path (report.dir, "supplementary-reanalysis-GSE9893.tex"))

## sets working directory back to project's root directory
setwd ("../..")

## * Quits
## sessionInfo ()
## q (save = "no")
