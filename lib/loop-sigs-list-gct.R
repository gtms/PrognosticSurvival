## loopSigsListGCT
## Gil Tomás, IRIBHM - ULB
## gil.tomas@ulb.ac.be

## this is an internal of testMSigDBgct

loopSigsListGCT <- function (gct,
                             surv,
                             sig.lst,
                             pca.method = "svd",
                             minGenes = 5,
                             doParallel = FALSE,
                             scale.matrix = FALSE) {
    ## only considers gene signatures with at least 'minGenes' present in the
    ## expression matrix
    isNA <- laply (sig.lst, function (sig) {
        length (intersect (sig, gct$row.descriptions)) < minGenes
    },
                   .parallel = doParallel)
    cox.dfr <- ldply (sig.lst[!isNA], function (sig) {  
        score <- binaryPC1Score (gct,
                                 sig,
                                 method = pca.method,
                                 scale.matrix = scale.matrix)
        data.frame (logRankTest (score$discrete,
                                 surv[!score$rm.samples]))
    },
                      .parallel = doParallel)
    names (cox.dfr)[1] <- "signature"
    cox.dfr
}
