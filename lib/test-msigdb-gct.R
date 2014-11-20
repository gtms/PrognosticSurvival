## testMSigDBgct

## Retrieves output of logrank tests for all MSigDB signatures, given a gct and
## an object of class survival.

testMSigDBgct <- function (gct,
                           msigdb.lst = NULL,
                           surv, # object of class Survival
                           minGenes = 5, # only tests for signatures that share
                                        # at least 'minGenes' genes with
                                        # expression profile
                           is.rnd = FALSE, # test with randomized MSigDB signatures?
                           pca.method = "svd",
                           scale.matrix = TRUE,
                           nCores = NULL) {
  if (is.null (nCores)) {
      nCores <- detectCores ()
  }
  if (nCores == 1) doParallel <- FALSE else doParallel <- TRUE
  registerDoMC (nCores)
  ## if msigdb.lst is NULL, takes loaded MSigDB C2 parsed signatures
  if (!is.null (msigdb.lst)) {
      rm (mSigDBc2)
      mSigDBc2 <- msigdb.lst
  }
  ## selects list of gene expression signatures
  if (is.rnd) {
      sig.lst <- mSigDBc2$rnd.sigs
  } else {
      sig.lst <- mSigDBc2$sigs
  }
  ## removes samples in gct for which no survival data is available
  toKeep <- !is.na (surv)
  gctRdx <- gct
  gctRdx$data <- gctRdx$data[, toKeep]
  survRdx <- surv
  survRdx <- surv[toKeep]
  ## performs log rank test for each signature in sig.lst
  loopSigsListGCT (gct = gctRdx,
                   surv = survRdx,
                   sig.lst = sig.lst,
                   pca.method = pca.method,
                   minGenes = minGenes,
                   doParallel = doParallel,
                   scale.matrix = scale.matrix)
}
