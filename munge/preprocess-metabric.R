## Preprocesses METABRIC data for dataset-specific analysis
## Vincent Detours, IRIBHM - ULB
## vdetours@ulb.ac.be

library (ProjectTemplate)
load.project ()

## * Merges 'discovery' and 'validation' metabric datasets
## ** loads data
dsc.dset <- loadDset ("metabric-discovery-set")
val.dset <- loadDset ("metabric-validation-set")

## ** merges data
metabric.dset <- list ()
## gct
metabric.dset$gct <- mergeGCTs (list (dsc.dset$gct, val.dset$gct))
## pheno
dsc.pheno <- dsc.dset$pheno
dsc.pheno$series <- "discovery"
val.pheno <- val.dset$pheno
val.pheno$series <- "validation"
metabric.dset$pheno <- rbind (dsc.pheno, val.pheno)
metabric.dset$pheno$os <- with (metabric.dset$pheno, Surv (os[, "time"], os[, "status"]))
metabric.dset$pheno$dss <- with (metabric.dset$pheno, Surv (dss[, "time"], dss[, "status"]))

## * Caches merged data
ProjectTemplate:::cache ("metabric.dset")
