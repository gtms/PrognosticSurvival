# Supplementary Code

This directory tree contains supplementary code supporting the analysis reported
in the article "Pervasive prognostic signals in the cancer transcriptome or why
association with outcome is not biologically informative".

Author: Gil Tomás <gil.tomas@ulb.ac.be>

URL: <https://owncloud.ulb.ac.be/public.php?service=files&t=506a90a16b44d0e8050f2e991b70c7d8>


## Contents
1. Preliminaries and license
2. Install
3. Run the analysis
4. Where things are


## 1. Preliminaries and license

The execution of this code requires a LINUX/UNIX environment with a working R
(version>=3.0) and TeX installations.  Its sole intent is to support the
findings reported in the quoted article.

This file is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This file is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this file.  If not, see <http://www.gnu.org/licenses/>.


## 2. Install

Prior requirements to the running of this software include a working R
(version>=3.0) and TeX installations.  Furthermore, R packages described in the
file `config/global.dcf` are also expected to be found in your system.  In
addition, the CRAN R package ProjectTemplate (version>=0.6) and
MicroarrayToolbox (available on <http://github.com/gtms/MicroarrayToolbox>) must
also be installed.

An R script located on `install/install-packages.R` can be executed to fill in these
requirements.  On a bash command line, enter:

```bash
R CMD BATCH install/install-packages.R
```


## 3. Run the analysis

This project runs within the R ProjectTemplate framework for automated data
analysis (<http://projecttemplate.net>).

Launch R at the root directory of the project, where this README.md file is
located, or set the working directory with the `setwd ()` command.

Then you need to run the following two lines of R code:

```R
library ("ProjectTemplate")
load.project ()
```

Once the second line of code is evaluated, a series of automated tasks will be
executed depending on the configurations declared in the `config/global.dcf`
file.  With the original configuration, these tasks include:
* Loading any R packages listed in the configuration file.
* Reading relevant datasets stored in `data` or `cache`.
* Pre-processing the data using the files in the `munge` directory.
* Executing the analysis of pre-processed data, yielding graphical output data.


## 4. Where things are

* Directories
    - _Configuration files_

      The analysis work-flow followed by ProjectTemplate is determined by the
      configuration flags found in the `config/global.dcf` file.  Depending on
      the `TRUE/FALSE` status of these flags, the function `load.project ()`
      may: load raw data into memory (flag `data_loading`); load pre-processed
      data into memory (flag `cache_loading`); pre-process raw data (flag
      `munging`); and load pre-determined libraries into memory (flag
      `load_libraries`).  For instance, once the raw data has been initially
      pre-processed and cached, you may find it desirable to turn the `munging`
      flag off and the `cache_loading` flag on.  This will allow for direct
      access to pre-processed data on your working environment once
      `load.project ()` is executed on later R sessions.
    - _Raw data_

      Raw data can be found in the `data` directory.  The `csv` directory
      contains the file `studies.csv`, which has information about all data-sets
      analyzed in this study.  The `rda` directory contains all data-sets stored
      on disk as `Rda` files.  The `sigs` directory contains biologically
      motivated gene expression signatures in `Rda` format, plus the 4722 MSigDB
      curated gene sets (collection v4.0, updated on May 31, 2013), as
      downloaded from <http://www.broadinstitute.org/gsea/msigdb>, in the `gmt`
      format.
    - _Pre-processed data_

      Pre-processed data, or cached data, can be found in the `cache` directory as
      `Rda` files.  These are the output of the processing of raw data with scripts
      located in the `munge` directory.
    - _Pre-processing scripts_

      Pre-processing scripts are located in the `munge` directory.
    - _Processing scripts_

      Processing scripts, carrying for the core of the analysis and producing
      its final graphical output, are located in the `src` directory.
    - _Supplementary Methods_

      Source `tex` files for the *Supplementary Methods* documentation can be
      found in the `reports` directory.  Upon processing with `pdflatex`, the
      corresponding `pdf` files will be located in the same directory.
    - The remaining directories should be self-explanatory.

* Output file types
    - `*.Rout`

      These are run logs, i.e. records of a particular computation run by a
      script.  These include non-graphical intermediate results, values of the
      random number generator seeds, and software package versions.
    - `*.pdf`

      These are graphical outputs.
    - `*.Rda`

      These are binary files storing pre-processed intermediate results that can
      be further scrutinized, should the user decide to fine-tune the analysis or
      push it in another direction.
