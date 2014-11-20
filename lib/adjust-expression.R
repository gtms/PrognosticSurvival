## adjustExpression

## Performs deconvolution of one (univariate) or more (multivariate) variables,
## defined in var.dfr, on the expression of a given gene, defined by expr.  This
## function is meant to be called by the global function adjustGCT.  The default
## method to be used for the resistant linear fit is "lqs", as defined in
## adjustGCT.  This function requires the package MASS to be loaded on the
## search path.

adjustExpression <- function (expr, # numeric vector of gene expression
                              var.dfr, # named data frame of variables to
                                       # deconvolve
                              method) {
    sampleNames <- names (expr)
    isNA <- is.na (expr) | Reduce ('|', lapply (var.dfr, is.na))
    frml <- reformulate (termlabels = names (var.dfr),
                         response = "expr")
    var.dfr$expr <- expr
    fit <- lqs (frml, method = method, data = var.dfr)
    res <- rep (NA, length (isNA))
    res[!isNA] <- fit$residuals
    out <- res + mean (expr, na.rm = TRUE)
    names (out) <- sampleNames
    return (out)
}
