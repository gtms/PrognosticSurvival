## Parses MSigDB C2 signatures
## Gil Tomás, IRIBHM - ULB
## gil.tomas@ulb.ac.be

## * Preamble
library (ProjectTemplate)
load.project ()

## * Code
## reads MSigDB c2 gene sets
lines <- readLines ("data/sigs/c2.all.v4.0.symbols.gmt")
sigs <- lapply (lines, function (l) do.call (c, strsplit (l, "\t")))
## names the sigs
sigs <- setNames (sigs, sapply (sigs, function (sig) sig[1]))
## discards the now spurious first two elements of each sig
sigs <- lapply (sigs, function (sig) sig[-c (1:2)])
## how many sigs?
length (sigs)
## have a look at the head of first 6 sigs
lapply (head (sigs), head)

## generates random signatures
set.seed (18042013)
all.genes <- unique (do.call (c, sigs))
system.time (rnd.sigs <- lapply (sigs, function (sig) {
    sample (all.genes, length (sig), replace = FALSE)
}))
## have a look at the head of the first 6 random sigs
lapply (head (rnd.sigs), head)

## packs, caches
mSigDBc2 <- list (sigs = sigs,
                  rnd.sigs = rnd.sigs)

ProjectTemplate:::cache ("mSigDBc2")

## * Visualization
## (uncomment to evaluate)
## size.distr <- sapply (sigs, length)
## summary (size.distr)
## hist (size.distr, breaks = 20)

## genes <- do.call (c, sigs)
## gene.count <- as.data.frame (table (genes),
##                              stringsAsFactors = FALSE)

## how many unique genes?
## dim (gene.count)
## gene.count <- gene.count[order (gene.count$Freq, decreasing = TRUE), ]
## gene.count[1:50, ]

## tail (gene.count, 50)

## how many genes are only represented once out of the 21046?
## sum (gene.count$Freq == 1)

## top.genes <- gene.count$genes[1:100]
## system.time (go <- getGO (top.genes))
## go.bp <- go[ontology == "BP", term]
## bp.count <- as.data.frame (table (go.bp),
##                            stringsAsFactors = FALSE)
## bp.count <- bp.count[order (bp.count$Freq, decreasing = TRUE), ]
## bp.count[1:50, ]

sessionInfo ()
