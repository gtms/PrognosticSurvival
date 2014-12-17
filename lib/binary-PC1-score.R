## binaryPC1Score
## defines a partition score method based on the first principal component
## Gil Tomás, IRIBHM - ULB
## gil.tomas@ulb.ac.be

## Given a gene expression matrix in the gct format (argument 'gct') and a gene
## signature (argument 'sig'), this function computes the first principal
## component of the signature across the patient's cohort.  It then returns a
## list containing: the loadings of each sample towards the first prinicpal
## component (object 'continuous'); a binary vector bisecting the samples
## according to the median value of the first principal component (object
## 'discrete'); a logical vector describing which samples were removed from the
## computation for not having enough genes profiled (object 'rm.samples').

binaryPC1Score <- function (gct,
                            sig,
                            method = "svdImpute",
                            best.pc = FALSE,
                            scale.matrix = FALSE) {
    removeNAcolumns <- function (gct,
                                 threshold = .1) { # minimum percentage of the
                                        # genes profiled in order to
                                        # retain the sample
        idd <- with (gct, apply (data, 2, function (x) sum (!is.na (x)) <= threshold * nrow (data)))
        list (gct = list (row.descriptions = gct$row.descriptions,
                  data = gct$data[, !idd, drop = FALSE]),
              rm.col = idd)
    }
    gct <- subsetGCTwithSig (gct = gct,
                             sig = sig)
    pruned.gct <- removeNAcolumns (gct, threshold = 0)
    gct <- pruned.gct$gct
    rm.samples <- pruned.gct$rm.col
    if (!is.matrix (gct$data)) {
        return ("None of the genes in the signature is profiled in the expression data.")
    }
    if (scale.matrix) scale <- "uv" else scale <- "none"
    res.pca <- pcaMethods:::pca (gct$data, # The package `amap' also contains a pca function
                                 method = method,
                                 scale = scale,
                                 verbose = FALSE)
    if (best.pc) {
        ## chooses PC best correlated with metagene among those capturing
        ## more than 5% of the variation in the dataset
        metagene <- computeSuperIndex (gct, sig)
        keep <- res.pca@R2 > .05
        cors <- abs (cor (loadings (res.pca), metagene))
        which.pc <- which.max (cors[keep])
        pc <- loadings (res.pca)[, which.pc]
        list (continuous = pc,
              discrete = as.numeric (pc > median (pc)),
              which.pc = which.pc,
              cor = max (cors),
              cor.test.pval = cor.test (loadings (res.pca)[, which.pc],
                  metagene)$p.value,
              rm.samples = rm.samples)
    } else {
        pc1 <- loadings (res.pca)[, 1]
        list (continuous = pc1,
              discrete = as.numeric (pc1 > median (pc1)),
              rm.samples = rm.samples)
    }
}
