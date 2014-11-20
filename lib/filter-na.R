## Given an object `obj' of class character or factor, replaces each of its
## instances matching `regexp' with NA, then returns `obj'.

filterNA <- function (obj,
                      regexp = "^[[:space:]]*(n(.?)a|non(.?)available)[[:space:]]*$") {
    if (!is.character (obj) & !is.factor (obj)) {
        stop ("filterNA: `obj' must be of class character or factor.")
    }
    if (is.null (obj)) return (NA)
    strg <- as.character (obj)
    strg[sapply (strg, nchar) == 0] <- NA
    strg[sapply (strg, grepl, pattern = regexp, ignore.case = TRUE)] <- NA
    if (is.factor (obj)) {
        out <- obj
        out[is.na (strg)] <- NA
    } else {
        out <- strg
    }
    return (out)
}

## test filterNA
## vec <- c ("N/A", " na", "Navailable", "National", "n+a",
## "non++available", "n+a.", ".n+a", "n/A ", " nyA", "nOn-avAilabLe ", "NONAVAILABLE")
## fac <- factor (vec)
## filterNA (vec)
## filterNA (fac)
