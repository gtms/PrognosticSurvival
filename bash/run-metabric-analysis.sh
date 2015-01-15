#!/bin/bash
# This script runs the src/metabric.R script

prjDir=~/ownCloud/projects/prognostic-survival

nice /usr/bin/R CMD BATCH ${prjDir}/src/metabric.R
echo "Finished running script src/metabric.R on $(date)."

exit
