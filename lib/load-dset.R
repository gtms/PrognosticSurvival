## loadDset
## Gil Tomás, IRIBHM - ULB
## gil.tomas@ulb.ac.be

## loads given study dataset into current environment

loadDset <- function (study) {
    studies.dfr <- read.csv2 ("data/csv/studies.csv",
                              stringsAsFactors = FALSE)
    studies.dtb <- as.data.table (studies.dfr)
    setkey (studies.dtb, "study")
    load (sprintf ("data/rda/%s", studies.dtb[study][, rda]))
    return (dset)
}
