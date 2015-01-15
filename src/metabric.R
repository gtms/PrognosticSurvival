## Down-sampling experiments with the METABRIC cohort
## Vincent Detours, IRIBHM - ULB
## vdetours@ulb.ac.be

## * Preamble
library (ProjectTemplate)
load.project ()

## registers parallel backend with the 'foreach' package
registerDoMC ()

## sets constants
rep <- 100 # number of replicate runs
event <- "dss" # endpoint considered in this analysis
set.seed (06081945) # sets the seed

## * Analysis
gct <- metabric.dset$gct
surv <- metabric.dset$pheno[[event]]

## ** survey sampling-related variability
nSamples <- 500 # number of samples to bootstrap each time

system.time (fracGenesSignifMetabric <- bootstrapMetabric (gct, nSamples, rep, surv))

## ** survey effect of cohort size
nSamples <- c (100, 200, 300, 500, 750, 1000, 1250, 1500, 1750)

fracGenesSignifMetabricCohortSize.lst <- lapply (nSamples, function (nb) {
    bootstrapMetabric (gct, nb, rep, surv)
})

names (fracGenesSignifMetabricCohortSize.lst) <- nSamples

## ** survey follow-up time
## define follow-up baseline follow-up times (in months)
nSamples <- 500 # number of samples to bootstrap each time
surv.mtx <- unclass (surv)
surv.time <- surv.mtx[, "time"]
surv.status <- surv.mtx[, "status"]
fut <- seq (25, 200, 25)
fut <- c (fut, max (surv.time, na.rm = TRUE))

system.time (fracGenesSignifMetabricFollowUp.lst <- lapply (fut, function (t) {
    surv.status[surv.time > t] <- 0 # events beyond horizon time not observed
    surv.time[surv.time > t] <- t # horizon time becomes ceiling time
    surv <- Surv (surv.time, surv.status)
    bootstrapMetabric (gct, nSamples, rep, surv)
}))

names (fracGenesSignifMetabricFollowUp.lst) <- fut

## ** survey ER status
nSamples <- 300
er <- metabric.dset$pheno$ER_IHC_status # We chose the IHC standard to call the ER status
fracERplus <- seq (0, 1, 0.1)

system.time (fracGenesSignifMetabricERplus.lst <- lapply (fracERplus, function (fep) {
    sapply (1:rep, function (n) {
        i <- c (sample (which (er == "pos"), nSamples * fep),
                sample (which (er == "neg"), nSamples * (1 - fep)))
        bootstrap.gct <- list ()
        bootstrap.gct$row.descriptions <- gct$row.descriptions
        bootstrap.gct$data <- gct$data[, i]
        bootstrap.surv <- surv[i]
        frac.signif.pval (bootstrap.gct, bootstrap.surv)
    })
}))

names (fracGenesSignifMetabricERplus.lst) <- fracERplus

## ** survey nodal status
nSamples <- 500
lnp <- as.numeric (metabric.dset$pheno$lymph_nodes_positive) > 0 # mumber of invaded lymph nodes
fracLNplus <- seq (0, 1, 0.1)

system.time (fracGenesSignifMetabricLNplus.lst <- lapply (fracLNplus, function (fnp) {
    sapply (1:rep, function (n) {
        i <- c (sample (which (lnp), nSamples * fnp),
                sample (which (!lnp), nSamples * (1 - fnp)))
        bootstrap.gct <- list ()
        bootstrap.gct$row.descriptions <- gct$row.descriptions
        bootstrap.gct$data <- gct$data[, i]
        bootstrap.surv <- surv[i]
        frac.signif.pval (bootstrap.gct, bootstrap.surv)
    })
}))

names (fracGenesSignifMetabricLNplus.lst) <- fracLNplus

## ** survey cellularity
nSamples <- 175
cll <- as.character (metabric.dset$pheno$cellularity)
cllLevels <- c ("low", "moderate", "high")

system.time (fracGenesSignifMetabricCellularity.lst <- lapply (cllLevels, function (lvl) {
    idx <- cll %in% lvl
    gct <- subsetGCTsamples (gct, idx)
    bootstrapMetabric (gct, nSamples, rep, surv)
}))

names (fracGenesSignifMetabricCellularity.lst) <- cllLevels

## * Packs, caches
metabricSims.lst <- list ("sampling.variability" = fracGenesSignifMetabric,
                          "cohort.size" = fracGenesSignifMetabricCohortSize.lst,
                          "follow.up.time" = fracGenesSignifMetabricFollowUp.lst,
                          "er.status" = fracGenesSignifMetabricERplus.lst,
                          "nodal.status" = fracGenesSignifMetabricLNplus.lst,
                          "cellularity" = fracGenesSignifMetabricCellularity.lst)

ProjectTemplate:::cache ("metabricSims.lst")

## * Quits
sessionInfo ()
q (save = "no")
