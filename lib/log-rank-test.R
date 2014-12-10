## logRankTest
## Gil Tomás, IRIBHM - ULB
## gil.tomas@ulb.ac.be

## fits a Cox proportional hazards regression model (logrank test), formulating
## survival (as given by 'surv.obj') as a function of 'score'

logRankTest <- function (score,
                         surv.obj, ...) {
  coxph.smm <- summary (coxph (surv.obj ~ score, ...))
  ci <- coxph.smm$conf.int
  if (coxph.smm$coef[1] < 0 & identical (range (score), c (0, 1))) {
    ci <- 1 / ci
    x <- ci[3]
    ci[3] <- ci[4]
    ci[4] <- x
    score = 1 - score ## not great formally for non-binary, non-normalized scores
    coef <- -coxph.smm$coef[1]
  }
  list (hazard.ratio = ci[1],
        ci.low = ci[3],
        ci.high = ci[4],
        n = coxph.smm$n,
        p.val = unname (coxph.smm$logtest[3]),
        ebeta = coxph.smm$coef[1],
        survival.index = log10 (unname (coxph.smm$logtest[3])))
}
