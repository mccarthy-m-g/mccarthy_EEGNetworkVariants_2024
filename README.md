# McCarthy MSc Thesis

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/495956659.svg)](https://zenodo.org/badge/latestdoi/495956659)

[![Synced with](https://img.shields.io/badge/Synced%20with-OSF-blue.svg)](https://osf.io/xztdk/)
<!-- badges: end -->

This repository contains the materials and source code for my MSc thesis, titled "Studying Network Variants With Electroencephalography". This project uses the [targets](https://books.ropensci.org/targets/) R package to maintain a reproducible, end-to-end workflow for the entire computational pipeline. The configuration and computational steps in this pipeline can be found in `_targets.R`.

## Repository Overview

- `manuscripts/` contains all source files and rendered documents made by our targets pipeline, including the thesis document and supplementary materials.

- `figures/` contains all figures made by our targets pipeline.

- `_targets.R` contains the [target script file](https://docs.ropensci.org/targets/reference/tar_script.html) to configure and define our targets pipeline. This is executed using the `run.R` or `run.sh` scripts.

- `_targets/` contains, locally, the data store of our targets pipeline; in addition to metadata on the status of each target, which is included online.

- `data/` contains, locally, the raw and cleaned EEG and functional connectivity data, participant descriptive data; in addition to the data file used to handle EEG preprocessing for each recording, which is included online.

- `R/` contains the custom user-defined functions used in our targets pipeline.

- `python/` contains the scripts we used interactively during the initial preprocessing stages of the project, which were later refactored into our targets pipeline using the [reticulate](https://rstudio.github.io/reticulate/) R package.

- `LICENSE.md` contains a copy of the licenses content and code licenses for this repository.

## License

### Content License

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC_BY--NC--ND_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

Except where otherwise noted, this project is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](https://creativecommons.org/licenses/by-nc-nd/4.0/). Please see `LICENSE.md` for a copy of this license.

You are free to:

**Share**---copy and redistribute the material in any medium or format.

Under the following terms:

**Attribution**---You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.

**NonCommercial**---You may not use the material for commercial purposes.

**NoDerivatives**---If you remix, transform, or build upon the material, you may not distribute the modified material.

### Code License

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

All source code in this project is licensed under the [MIT License](https://opensource.org/licenses/MIT). Please see `LICENSE.md` for a copy of this license.

You are free to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so.

Under the following terms:

**Attribution**---You must include the following copyright notice "Copyright (c) 2022 Michael McCarthy" and provide a link to the license in all copies or substantial portions of the Software.
