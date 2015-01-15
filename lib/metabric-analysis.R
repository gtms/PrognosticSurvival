## Functions supporting the METABRIC dataset analysis
## Vincent Detours, IRIBHM - ULB
## vdetours@ulb.ac.be

## * Functions to compute fraction of genes associated with outcome
## ** get.cox.lr.pval
## given matrix of expression (mtx) and survival object (s), returns vector
## of p-value of individual gene associations with outcome
get.cox.lr.pval <- function (mtx, s) {
    aaply (mtx, 1,
           function (x) summary (coxph(s ~ x))$logtest[['pvalue']],
           .parallel = TRUE)
}

## ** frac.signif.pval
## given gct and survival object (s), returns fraction of genes associated with
## outcome
frac.signif.pval <- function (gct, s, p = 0.05) {
    mtx <- gct$data
    pval <- get.cox.lr.pval (mtx, s)
    sum (pval < p) / length (pval)
}

## ** bootstrapMetabric
## given a gct, a number of samples (nSamples) and a number or repetitions
## (rep), returns a vector with the fraction of genes associated with outcome
## within 'rep' bootstrappings of 'nSamples' withing the gct
bootstrapMetabric <- function (gct, nSamples, rep, surv) {
    sapply (1:rep, function (n) {
        i <- sample (1:ncol (gct$data), nSamples)
        bootstrap.gct <- list ()
        bootstrap.gct$row.descriptions <- gct$row.descriptions
        bootstrap.gct$data <- gct$data[, i]
        bootstrap.surv <- surv[i]
        frac.signif.pval (bootstrap.gct, bootstrap.surv)
    })
}

## * Functions to display analysis results
## ** transparent.rgb
transparent.rgb <- function(col, alpha = 30) {
    tmp <- c (col2rgb (col), alpha, 255)
    names (tmp) <- c ("red", "green", "blue", "alpha", "maxColorValue")
    do.call ("rgb", as.list (tmp))
}

## ** myboxplot
myboxplot <- function  (data,
                        yat = NULL,
                        xat = NULL,
                        xlab = "",
                        ylab = "",
                        main = "",
                        file = "",
                        col = "black", ...) {
    if (file != "") {
        do.call (gsub (".*([a-z]+{3})$", "\\1", file), list (file))
        par (mar = c (6.1, 6.5, 4.1, 1.1))
    }

    boxplot (data,
             outline = FALSE,
             frame = TRUE,
             lwd = 1,
             cex.axis = 2.5,
             cex.lab = 2.5, ...)
    stripchart (as.data.frame (data),
                add = TRUE,
                method = "jitter",
                pch = 19,
                vertical = TRUE,
                cex = 2,
                col = transparent.rgb (col))

    mtext (side = 1, text = xlab, line = 4, cex = 2.5, padj = 0.25)
    mtext (side = 2, text = ylab, line = 4, cex = 2.5)

    if (file != "") {
        dev.off()
    }
}
