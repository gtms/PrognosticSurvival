## adjustGCT
## Gil Tomás, IRIBHM - ULB
## gil.tomas@ulb.ac.be

## Performs deconvolution of one (univariate) or more (multivariate) variables,
## defined in var.dfr, on a per-gene basis across the expression matrix of a
## given gct object.  The default method to be used for the resistant linear fit
## is "lqs".  This function requires the packages MASS, plyr and doMC to be
## installed.

adjustGCT <- function (gct,
                       var.dfr, # named data frame of variables to deconvolve
                       method = "lqs",
                       is.control = FALSE, # whether to shuffle the variable
                                        # indices or not
                       nCores = NULL) {
    ## determines number of machine cores to be employed during computation
    if (is.null (nCores)) {
        nCores <- detectCores ()
    }
    if (nCores == 1) doParallel <- FALSE else doParallel <- TRUE
    registerDoMC (nCores)
    ## whether to shuffle variable indices of not
    if (is.control) var.dfr <- as.data.frame (sapply (var.dfr, sample))
    mtx <- gct$data
    adj.mtx <- aaply (mtx, 1, adjustExpression,
                      var.dfr = var.dfr,
                      method = method,
                      .parallel = doParallel)
    gct$data <- adj.mtx
    return (gct)
}
