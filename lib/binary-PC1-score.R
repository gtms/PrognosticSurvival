## binaryPC1Score
## defines a partition score method based on the first principal component

binaryPC1Score <- function (gct,
                            sig,
                            method = "svd",
                            scale.matrix = FALSE,
                            best.pc = FALSE,
                            removeNAcolumns = FALSE) {
    ## subsets gct with signature
    gct <- subsetGCTwithSig (gct = gct,
                             sig = sig)
    if (removeNAcolumns) {
        pruned.gct <- removeNAcolumns (gct, threshold = 0)
        gct <- pruned.gct$gct
        rm.samples <- pruned.gct$rm.col
    } else {
        rm.samples <- rep (NA, ncol (gct$data))
    }        
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
