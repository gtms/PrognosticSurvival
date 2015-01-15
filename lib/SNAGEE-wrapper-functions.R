## Wrapper functions for the qualStudy function of the SNAGEE package.  This
## function returns a signal/noise quality metric based on gene-gene
## correlations across several expression profiles.
## Gil Tomás, IRIBHM - ULB
## gil.tomas@ulb.ac.be

## computeQualStudy
computeQualStudy <- function (study) {
    dset <- loadDset (study)
    d <- gct2qualStudy (dset$gct)
    qualStudy (d)
}

## gct2qualStudy
## transforms gct into input for qualStudy
gct2qualStudy <- function (gct) {
    rd <- as.character (gct$row.descriptions)
    gene.symb <- ifelse (nchar (rd) == 0, "none", rd)
    entrez.ids <- laply (mget (gene.symb, org.Hs.egSYMBOL2EG, ifnotfound = NA),
                         function (x) as.numeric (x[1]))
    return (list (genes = entrez.ids,
                  data = gct$data))
}
