## Computes fraction of MSigDB C2 signatures associated with outcome
## across all studies.

## * Preamble
library (ProjectTemplate)
load.project ()

## registers parallel backend with the 'foreach' package
registerDoMC ()

## * Defines script-specific functions
## ** extractFracSignif
extractFracSignif <- function (dfr, signifCutoff = .05) {
    with (dfr, sum (p.val < signifCutoff, na.rm = TRUE) / sum (!is.na (p.val)))
}

## ** testSingleGene
testSingleGene <- function (sc,
                            surv.obj) {
    ## returns result of log rank test
    data.frame (logRankTest (score = sc,
                             surv = surv.obj))
}

## ** testAllGenes
testAllGenes <- function (gct,
                          surv.obj) {
    ## retains only individuals with non-absent survival data
    keep <- !is.na (surv.obj)
    mtx <- gct$data[, keep]
    surv.obj <- surv.obj[keep]
    cox.dfr <- adply (mtx, 1, testSingleGene,
                      surv = surv.obj,
                      .parallel = TRUE)
    ## computes q-values
    ## cox.dfr$q.val <- qvalue (cox.dfr$p.val)$qvalues
    cox.dfr$q.val <- p.adjust (cox.dfr$p.val, method = "fdr")
    cox.dfr[, 1] <- as.character (cox.dfr[, 1])
    names (cox.dfr)[1] <- "ID"
    ## returns data frame
    cbind (gene = gct$row.descriptions,
           cox.dfr)
}

## ** computeFracGenesSignif
computeFracGenesSignif <- function (study,
                                    restrict.event = FALSE,
                                    nCores = NULL,
                                    print.token = FALSE) {
    ## if nCores is not user specified, use all available cores in the machine
    if (is.null (nCores)) {
        nCores <- detectCores ()
    }
    ## entry token
    if (print.token) {
        message (sprintf ("Now doing %s . . . ", study), appendLF = FALSE)
    }
    ## loads study
    dset <- loadDset (study)
    ## if restrict.event, select which event to use
    available.events <- names (dset$pheno)[sapply (dset$pheno, class) == "Surv"]
    if (restrict.event) {
        all.events <- c ("dss", "os", "dmfs", "dfs")
        events <- all.events[all.events %in% available.events][1]
    } else {
        events <- available.events
    }
    ## loops over list of events
    allOut.lst <- lapply (events, function (event) {
        surv = dset$pheno[[event]]
        ## tests association of single genes with outcome
        res.dfr <- testAllGenes (gct = dset$gct,
                                 surv.obj = surv)
        res.dfr$event <- event
        res.dfr
    })
    ## packs log rank tests and fraction of significant tests per event
    logRankTests.dtb <- data.table (do.call (rbind, allOut.lst))
    logRankTests.dtb[, study := study]
    fracGenesSignif.dtb = logRankTests.dtb[,
        list (fracGenesSignif = extractFracSignif (.SD)),
        by = event]
    fracGenesSignif.dtb[, study := study]
    ## exit token
    if (print.token) {
        message (sprintf ("%s done!\n", study))
    }
    ## return
    list (logRanktests = logRankTests.dtb,
          fracGenesSignif = fracGenesSignif.dtb)
}


## * Reads in data
## ** reads in studies.csv2
studies.dfr <- read.csv2 ("data/csv/studies.csv",
                          stringsAsFactors = FALSE)

## * Performs computation
## ** computes fraction of genes significantly associated with outcome across all studies
## loops computeFracGenesSignif over all studies
## COMPUTATIONALLY INTENSIVE!
## With 16 cores:
##      user    system   elapsed 
## 21391.551  1311.320  5039.397 (roughly one hour and 25 minutes)
system.time (fracGenesSignif.lst <- setNames (llply (studies.dfr$study,
                                                     computeFracGenesSignif,
                                                     print.token = TRUE),
                                              studies.dfr$study))

fracGenesSignif.dtb <- do.call (rbind, lapply (fracGenesSignif.lst, function (lst) lst$fracGenesSignif))
fracGenesSignif.dfr <- as.data.frame (fracGenesSignif.dtb)

## * Caches results
ProjectTemplate:::cache ("fracGenesSignif.lst")
ProjectTemplate:::cache ("fracGenesSignif.dfr")
