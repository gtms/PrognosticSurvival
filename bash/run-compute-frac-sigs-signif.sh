#!/bin/bash
# This script runs the src/compute-frac-sigs-signif.R script

prjDir=~/ownCloud/projects/prognostic-survival

nice /usr/bin/R CMD BATCH ${prjDir}/munge/compute-frac-sigs-signif.R
echo "Done."

exit
