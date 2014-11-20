## computeSuperIndex
## computes super index from gct, given gene signature 

computeSuperIndex <- function (gct,
                               sig, # character vector of gene symbols
                               is.probeset = FALSE) {
    ifelse (is.probeset,
            .id <- rownames (gct$data) %in% sig,
            .id <- gct$row.descriptions %in% sig)
    if (sum (.id) == 0) return (NA)
    if (sum (.id) == 1) {
        gct$data[.id, ]
    } else {
        mtx <- gct$data[.id, , drop = FALSE]
        ## Should there be duplicated row names, collapse them by median
        ## expression
        if (sum (duplicated (gct$row.descriptions[.id])) != 0) {
            mtx <- as.matrix (apply (mtx, 2, function (sample) {
                tapply (sample, gct$row.descriptions[.id], median, na.rm = TRUE)
            }))
        }
        mtx.median <- median (mtx, na.rm = TRUE)
        if (mtx.median == 0) mtx.median <- 1
        mtx <- t (apply (mtx, 1, function (x) {
            x * mtx.median / ifelse (median (x, na.rm = TRUE) == 0, 1, median (x, na.rm = TRUE))
        }))
        apply (mtx, 2, function (y) median (y, na.rm = TRUE))
    }
}
