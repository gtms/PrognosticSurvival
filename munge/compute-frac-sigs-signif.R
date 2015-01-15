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

## ** computeFracSignif
computeFracSignif <- function (study,
                               restrict.event = FALSE,
                               discrete.pc1 = FALSE,
                               nCores = NULL,
                               return.all.tests = FALSE,
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
        ## defines res.lst
        res.lst <- list ()
        ## computes logRankTests for all MSigDB C2 signatures
        res.lst$msigdb <- testMSigDBgct (gct = dset$gct,
                                         surv = surv,
                                         discrete.pc1 = discrete.pc1,
                                         scale.matrix = TRUE,
                                         nCores = nCores)
        res.lst$random <- testMSigDBgct (gct = dset$gct,
                                         surv = surv,
                                         discrete.pc1 = discrete.pc1,
                                         is.rnd = TRUE,
                                         scale.matrix = TRUE,
                                         nCores = nCores)
        if (!return.all.tests) {
            out.lst <- lapply (res.lst, extractFracSignif)
        } else {
            out.lst <- res.lst
        }
        out.lst$event <- event
        out.lst$study <- study
        out.lst
    })
    ## exit token
    if (print.token) {
        message (sprintf ("%s done!\n", study))
    }
    ## returns list or data frame, depending on return.all.tests
    if (!return.all.tests) {
        Reduce (rbind, lapply (allOut.lst, as.data.frame))
    } else {
        allOut.lst
    }
}

## * Reads in data
## ** reads in studies.csv2
studies.dfr <- read.csv2 ("data/csv/studies.csv")

## * Performs computation
## ** computes fraction of significant tests across all studies
## loops computeFracSignif over all studies
## COMPUTATIONALLY INTENSIVE!
## With 16 cores:
##      user    system   elapsed 
## 150432.18  10046.56  13398.25 (almost 4 hours)
## system.time (fracSignif.dfr <- ldply (studies.dfr$study,
##                                       computeFracSignif,
##                                       print.token = TRUE))

system.time (fracSignifDiscretePC1.dfr <- ldply (studies.dfr$study,
                                                 computeFracSignif,
                                                 discrete.pc1 = TRUE,
                                                 print.token = TRUE))

## system.time (fracSignif.lst <- llply (studies.dfr$study,
##                                       computeFracSignif,
##                                       return.all.tests = TRUE,
##                                       print.token = TRUE))

## * Caches results
## ProjectTemplate:::cache ("fracSignif.lst")
## ProjectTemplate:::cache ("fracSignif.dfr")
ProjectTemplate:::cache ("fracSignifDiscretePC1.dfr")
