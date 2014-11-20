## removeNAcolumns
## Given a gct, returns the gct only with the samples successfully profiled for a
## minimum percentage of 'threshold' genes.

removeNAcolumns <- function (gct,
                             threshold = .9) { # minimum percentage of the
                                        # genes profiled in order to
                                        # retain the sample
    idd <- with (gct, apply (data, 2, function (x) sum (!is.na (x)) <= threshold * nrow (data)))
    list (gct = list (row.descriptions = gct$row.descriptions,
              data = gct$data[, !idd, drop = FALSE]),
          rm.col = idd)
}
