## parsePheno

## Bioconductor expression sets (eSet) may harbour a phenoData object
## that contains information on experimental phenotypes recorded.  This object
## can be accessed with the 'pData' method and yields a data frame containing a
## number of clinical variables of interest in addition to other clinically
## irrelevant information, coded as factors.  The 'parsePheno' function takes
## this data frame as input and returns a cleaned data table, with the
## clinically relevant variables adequately encoded as objects of class 'factor'
## or 'numeric'.

parsePheno <- function (dfr,
                        sep = ": ") {
    require (data.table)
    trim.dfr <- dfr[, grepl ("characteristics", colnames (dfr)), drop = FALSE]
    trim.dfr <- as.data.frame (lapply (trim.dfr, as.character))
    ## extract colnames
    extractClnms <- function (vec, sep) {
        splt.lst <- lapply (vec, strsplit, split = sep)
        firstSplt.lst <- lapply (splt.lst, function (y) {
            if (length (y[[1]]) != 0) y[[1]][[1]][1] else NULL
        })
        unique (unlist (firstSplt.lst))
    }
    ## fetches list of variables coaxed in each column
    clnms.lst <- lapply (trim.dfr, extractClnms, sep = sep)
    ## determines set of unique clinical variables in the data frame
    clnms <- unique (unlist (clnms.lst))
    ## fills in variables
    seekVar <- function (vec, clnms, sep) {
        vec <- unlist (vec)
        clnms <- sapply (clnms, escapeMeta)
        var.lst <- lapply (clnms, function (strg) {
            logVec <- grepl (sprintf ("^%s%s", strg, sep), vec)
            if (sum (logVec) == 0) {
                NA
            } else {
                strsplit (vec[logVec], sep)[[1]][2]
            }
        })
        data.frame (setNames (var.lst, clnms))
    }
    ## loops over seekVar to build data frame of clinical variables
    clin.dfr <- adply (trim.dfr, 1, seekVar, clnms = clnms, sep = sep)
    ## clin.dfr <- clin.dfr[, gsub ("[[:space:]]", ".", clnms)]
    clin.dfr <- clin.dfr[, !grepl ("^characteristics", names (clin.dfr))]
    ## converts variables to factors
    clin.dfr <- data.frame (lapply (clin.dfr, as.factor))
    ## filters NA variables
    clin.dfr <- data.frame (lapply (clin.dfr, filterNA))
    ## transforms potential numeric variables into numeric type
    clin.dfr <- transformDataFrameVars (clin.dfr)
    ## adds sample names
    rownames (clin.dfr) <- rownames (dfr)
    ## returns data table
    data.table (clin.dfr, keep.rownames = TRUE)
}
