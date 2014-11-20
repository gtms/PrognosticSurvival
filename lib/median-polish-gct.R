## medianPolish
## performs median polish of a gct

medianPolish <- function (gct) {
    if (sum (is.na (gct$data) != 0)) {
        gct <- pruneFeatures (gct) # MicroarrayToolbox
        gct$data <- impute.knn (gct$data)$data
    }
    gct$data <- medpolish (gct$data)$residuals
    gct
}
