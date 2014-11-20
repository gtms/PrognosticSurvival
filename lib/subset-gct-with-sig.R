## subsetGCTwithSig
## subsets gct data with signatures's genes

subsetGCTwithSig <- function (gct,
                              sig) { # character vector of gene symbols
    idd <- gct$row.descriptions %in% sig
    list (row.descriptions = gct$row.descriptions[idd],
          data = gct$data[idd, , drop = FALSE])
}
