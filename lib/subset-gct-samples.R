## subsetGCTsamples
## Gil Tomás, IRIBHM - ULB
## gil.tomas@ulb.ac.be

## subsets gct's samples with a logical vector

subsetGCTsamples <- function (gct,
                              log.vec) { # logical vector with length
                                        # equal to the number of
                                        # samples in the gct
    list (row.descriptions = gct$row.descriptions,
          data = gct$data[, log.vec])
}
