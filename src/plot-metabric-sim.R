## Produces figure 2 of the paper, plotting the output of the script
## munge/metabric.R
## output saved in graphs/metabric-sims.R
## Gil Tomás and Vincent Detours, IRIBHM - ULB
## gil.tomas@ulb.ac.be

## * Preamble
library (ProjectTemplate)
load.project ()

## defines script-specific functions
transparent.rgb <- function (col, alpha = 30) {
    tmp <- c (col2rgb (col), alpha, 255)
    names (tmp) <- c ("red", "green", "blue", "alpha", "maxColorValue")
    do.call ("rgb", as.list (tmp))
}

myBoxplot <- function (data,
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
    par (mar = c (6.1, 6.5, 4.1, 1.1))
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
                ## cex = 2,
                col = transparent.rgb (col))
    mtext (side = 1,
           text = xlab,
           line = 4,
           cex = 2.5,
           padj = 0.25)
    mtext (side = 2,
           text = ylab,
           line = 4,
           cex = 2.5)
    if (file != "") {
        dev.off ()
    }
}

## defines script-specific variables
frac.label <- "% outcome-associated genes"

## * Preprocess
## multiply all values by 100 to yield percentages
plt.lst <- setNames (lapply (metabricSims.lst, function (lst) {
    setNames (lapply (lst, function (vec) vec * 100), names (lst))
}),
                     names (metabricSims.lst))

## unlist plt.lst[["sampling.variability"]]
plt.lst[["sampling.variability"]] <- unlist (plt.lst[["sampling.variability"]])

## * Plot
## defines layout of device; initializes device
svg (file = "graphs/metabric-sims.svg")
par (mar = c (6.1, 6.5, 4.1, 1.1))

## sets layout (2 x 3)
layout (mat = matrix (1:6, ncol = 3, byrow = TRUE))

## panel A
## sampling.variability
hist (plt.lst[["sampling.variability"]],
      xlab = frac.label,
      main = "",
cex.lab = 2.5,
cex.axis = 2.5)

## panel B
## cohort.size
myBoxplot (plt.lst[["cohort.size"]],
           xlab = "Number of patients",
           ylab = frac.label,
           xaxt = "n")
## file = "boxplot-metabric-N.svg")
## plot sample size (on top, first 4)
sample.size.top <- as.numeric (names (plt.lst[["cohort.size"]][1:4]))
text (x = 1:4 - .15,
      y = max (do.call (c, plt.lst[["cohort.size"]])) - 2,
      labels = sprintf ("n = %d", sample.size.top),
      srt = 90, # 90 degrees
      pos = 1) # right justififed
sample.size.bot <- as.numeric (names (plt.lst[["cohort.size"]][5:9]))
text (x = 5:9 - .15,
      y = min (do.call (c, plt.lst[["cohort.size"]])) + 1,
      labels = sprintf ("n = %d", sample.size.bot),
      srt = 90, # 90 degrees
      pos = 4) # left justififed

## panel C
## follow.up.time
myBoxplot (plt.lst[["follow.up.time"]],
           xlab = "Maximum follow-up (years)",
           ylab = frac.label,
           xaxt = "n")
## file = "boxplot-metabric-follow-up.svg")
fut.top <- as.numeric (names (plt.lst[["follow.up.time"]][1])) / 12
text (x = 1 - .15,
      y = max (do.call (c, plt.lst[["follow.up.time"]])) - 2,
      labels = sprintf ("%.1f yrs", fut.top),
      srt = 90, # 90 degrees
      pos = 1) # right justififed
fut.bot <- as.numeric (names (plt.lst[["follow.up.time"]][2:9])) / 12
text (x = 2:9 - .15,
      y = min (do.call (c, plt.lst[["follow.up.time"]])),
      labels = sprintf ("%.1f yrs", fut.bot),
      srt = 90, # 90 degrees
      pos = 4) # left justififed

## panel D
## er.status
myBoxplot (plt.lst[["er.status"]],
           xlab = "% ER+ tumors",
           ylab = frac.label)
## file = "boxplot-metabric-er.svg")

## panel E
## nodal.status
myBoxplot (plt.lst[["nodal.status"]],
           xlab = "% node positive patients",
           ylab = frac.label)
## file = "boxplot-metabric-node.svg")

## panel F
## cellularity
## myBoxplot (plt.lst[["cellularity"]],
##            xlab="Cellularity",
##            ylab=frac.label)
## file="boxplot-metabric-cellularity-class.svg")

## closes device
dev.off ()
