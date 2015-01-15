## Reanalysis of GSE9893: quantifying the effect of a noramlization artifact in
## prognostic content
## Gil Tomás & Vincent Detours, IRIBHM - ULB
## vdetours@ulb.ac.be

## * Preamble
library (ProjectTemplate)
load.project ()

## setup location for graphical output
graph.path <- "graphs/reanalysis-GSE9893"
dir.create (graph.path)

## loads study GSE9893-breast
GSE9893.dset <- loadDset ("GSE9893-breast")

## loads studies.dfr
studies.dfr <- read.csv2 ("data/csv/studies.csv", stringsAsFactors = FALSE)

## registers parallel backend with the 'foreach' package
registerDoMC ()

## * Analysis
## ** computes signal/noise quality metric for all datasets for comparison
system.time (studiesSNQual <- laply (studies.dfr$study, computeQualStudy, .parallel = TRUE))

names (studiesSNQual) <- studies.dfr$study

## ** demonstration of the batch effect and of its consequences
## reads and cleans GSE9893-gpr.csv
gpr.dfr <- read.csv2 ("data/raw-data-reanalysis-GSE9893/GSE9893-gpr.csv")
colnames (gpr.dfr) <- sapply (gpr.dfr[1, ], gsub,
                              pattern = "^(.*)=(.*)",
                              replacement = "\\1")
colnames (gpr.dfr)[1] <- "gsm"
gpr.dfr$gsm <- gsub (".gpr", "", gpr.dfr$gsm)
## removes constant variables
gpr.dfr <- gpr.dfr[sapply (gpr.dfr, function (x) length (unique (x)) > 1)]
## to data table
gpr.dtb <- data.table (gpr.dfr, key = "gsm")
## reorders according to GSE9893.dset$pheno
gpr.dtb <- gpr.dtb[rownames (GSE9893.dset$pheno)]
## whenever appropriate, removes everything before '=' (included) from each column
gpr.dtb <- gpr.dtb[, lapply (.SD, gsub, pattern = "^(.*)=(.*)", replacement = "\\2")]
## transform variables to numeric
numVar <- c ("Temperature", "PMTGain", "ScanPower", "LaserPower")
for (j in numVar) set (gpr.dtb, j = j, value = as.numeric (gpr.dtb[[j]]))
## creates new variable defined by the presence of the string 'bordeaux' in the
## 'Settings' variable
gpr.dtb[, Bordeaux := grepl ("bordeaux", Settings)]
## creates new variable by extracting the series number referenced in the
## variable 'ImageFiles'
gpr.dtb[, Series := gsub (".*rie[[:space:]]*([0-9]+).*", "\\1", ImageFiles)]
gpr.dtb$Series <- as.numeric (gpr.dtb$Series)

## arrays were scanned within two discrete time intervals
gpr.dtb[, Date := as.Date (DateTime, "%Y/%m/%d")]
gpr.dtb[, Time := as.numeric (Date)]
gpr.dtb[, Year := as.numeric (gsub ("^([0-9]+)-(.*)", "\\1", Date))]

svg (file.path (graph.path, "hist-time.svg"))
par (mar = c (6.1, 6.5, 4.1, 1.1))
with (gpr.dtb, hist (Time - min (Time),
                     breaks = 50,
                     main = "",
                     xlab = "Scan date [days]",
                     cex.lab = 2.5,
                     cex.axis = 2.5))
dev.off ()

## scan data is associated with outcome
summary (coxph (GSE9893.dset$pheno$os ~ gpr.dtb$Time))
summary (coxph (GSE9893.dset$pheno$os ~ gpr.dtb$Year))

svg (file.path (graph.path, "km-time.svg"))
par (mar = c (6.1, 6.5, 4.1, 1.1))
plot (survfit (GSE9893.dset$pheno$os ~ gpr.dtb$Year),
      xlab = "Time [months]",
      ylab = "Survival",
      cex.lab = 2.5,
      cex.axis = 2.5,
      col = 1:2)
dev.off ()

## there is a massive batch effect associated with scan year
pcaRes <- prcomp (t (GSE9893.dset$gct$data))
svg (file.path (graph.path, "PCA-loadings.svg"))
par (mar = c (6.1, 6.5, 4.1, 1.1))
plot (pcaRes,
      cex.lab = 2.5,
      cex.axis = 2.5)
dev.off ()

pt <- rep (1, nrow (GSE9893.dset$gct$data))
pt[GSE9893.dset$pheno$status == "RF"] <- 3
col <- rep ("black", nrow (GSE9893.dset$gct$data))
col[gpr.dtb$Year == "2006"] <- "red"
svg (file.path (graph.path, "PCA.svg"))
par (mar = c (6.1, 6.5, 4.1, 1.1))
plot (pcaRes$x[, 1], pcaRes$x[, 2], col = col, pch = pt,
      xlab = paste ("PC1 (", round (100 * pcaRes$sdev[1] ^ 2 / sum (pcaRes$sdev ^ 2),
          d = 1), "%)", sep = ""),
      ylab = paste ("PC2 (", round (100 * pcaRes$sdev[2] ^ 2 / sum (pcaRes$sdev ^ 2),
          d = 1), "%)", sep = ""),
      cex.lab = 2.5,
      cex.axis = 2.5,
      cex = 2)
dev.off()

## the batch effect is also visible from array-wise expression distributions
svg (file.path (graph.path, "boxplot-expression.svg"),
     height = 6,
     width = 13)
par (mar = c (6.1, 6.5, 4.1, 1.1))
boxplot (GSE9893.dset$gct$data,
         range = 0,
         col = col,
         xlab = "Samples",
         ylab = "Expression",
         pch = ".",
         lwd = 0.5)
dev.off ()

## ** renormalization from raw data
## signal/noise quality, but decreases prognostic fraction

## reads gpr files (these are one-color arrays scans)
gpr.dir <- "data/raw-data-reanalysis-GSE9893/GSE9893-gpr"
gpr.files <- file.path (gpr.dir, list.files (gpr.dir))
    
ff <- function (x) as.numeric(x$Flags > -99)
Cy5 <- "F635 Mean" # fools limma to read single clolor genepix file
RG <- read.maimages (gpr.files, source = "genepix.median",
                     wt.fun = ff,
                     columns = list (R = Cy5, G = Cy5))
RG$G <- NULL

## quantile-nomalize arrays
norm <- normalizeBetweenArrays (RG$R, method = "quantile")
norm <- log2 (norm)

## reorder columns
colnames (norm) <- gsub ("(.*)/(.*)$", "\\2", colnames (norm))
norm <- norm[, colnames (GSE9893.dset$gct$data)]

## gene-wise average
norm <- aggregate (norm, by = list (RG$genes$Name), mean)
rownames (norm) <- norm[, 1]
norm <- as.matrix (norm[, -1])

## PCA: no more batch effect
pcaResNorm <- prcomp (t (norm))
svg (file.path (graph.path, "PCA-loadings-norm.svg"))
par (mar = c (6.1, 6.5, 4.1, 1.1))
plot (pcaResNorm,
      cex.lab = 2.5,
      cex.axis = 2.5,
      cex = 2.0)
dev.off ()

svg (file.path (graph.path, "PCA-norm.svg"))
par (mar = c (6.1, 6.5, 4.1, 1.1))
plot (pcaResNorm$x[, 1],
      pcaResNorm$x[,2],
      col = col, pch = pt,
      xlab = paste ("PC1 (", round (100 * pcaResNorm$sdev[1] ^ 2 / sum (pcaResNorm$sdev ^ 2),
          d = 1), "%)", sep = ""),
      ylab = paste ("PC2 (", round (100 * pcaResNorm$sdev[2] ^ 2 / sum (pcaResNorm$sdev ^ 2),
          d = 1), "%)", sep = ""),
      cex.lab = 2.5,
      cex.axis = 2.5,
      cex = 2)
dev.off ()

## obviously, array expression distribution are now indentical
svg (file.path (graph.path, "boxplot-norm-expression.svg"),
     height = 6,
     width = 13)
par (mar = c (6.1, 6.5, 4.1, 1.1))
boxplot (norm,
         range = 0,
         col = col,
         xlab = "Samples",
         ylab = "Expression",
         pch = ".",
         lwd = 0.5)
dev.off ()

## concomitantly, prognostic fraction is sharpely reduced:
frac.signif.pval (gct = list (data = norm),
                  GSE9893.dset$pheno$os)


## yet signal/noise quality metric is improved and is now within the IQR of the
## surveyed datasets
norm.gct <- list (row.descriptions = rownames (norm),
                  data = norm)
(GSE9893NormSNQual <- qualStudy (gct2qualStudy (norm.gct)))

quantile (studiesSNQual, c(0.25, 0.75))

## displays distribution of signal/noise quality metrics
svg (file.path (graph.path, "hist-SN-metrics.svg"))
par (mar = c (6.1, 6.5, 4.1, 1.1))
hist (studiesSNQual,
      breaks = 20,
      main = "",
      xlab = "S/N metrics",
      cex.lab = 2.5,
      cex.axis = 2.5)
abline (v = studiesSNQual[["GSE9893-breast"]],
        col = "red",
        lwd = 1.5)
abline (v = GSE9893NormSNQual,
        col = "green",
        lwd = 1.5)
dev.off ()

## * Quits
sessionInfo ()
q (save = "no")
