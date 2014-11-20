## Takes data frame as input and, wherever sensible, transforms factor or
## character variables into numeric type.  Returns the transformed data frame.

transformDataFrameVars <- function (dfr) {
    ## fac2num: converts factor to numeric
    fac2num <- function (fac) {
        as.numeric (levels (fac))[as.integer (fac)]
    }
    ## which variables can potentially be converted to class numeric?
    toNum <- vapply (dfr, function (x) {
        sum (grepl ("[[:alpha:]]", x), na.rm = TRUE) == 0
    },
                     FUN.VALUE = TRUE)
    ## performs conversion
    dfr[toNum] <- sapply (dfr[toNum], function (x) {
        if (class (x) == "factor") {
            fac2num (x)
        } else {
            if (class (x) == "character") as.numeric (x)
        }
    })
    ## returns data frame
    dfr
}
