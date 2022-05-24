# McCarthy MSc Thesis

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/495956659.svg)](https://zenodo.org/badge/latestdoi/495956659)

[![Synced with](https://img.shields.io/badge/Synced%20with-OSF-blue.svg)](https://osf.io/xztdk/)
<!-- badges: end -->

<!--
## TODO:

-   set up renv package and docker
-   consider adding abstract
-->

## About

This repository contains the materials and source code for my MSc thesis, tentatively titled "Studying network variants with electroencephalography".

## Reproducibility

This project comes with a reproducible environment that includes everything needed to run its code. The reproducible environment is built using [Docker](https://www.docker.com), which packages all of the software dependencies of the project (i.e., the operating system, system libraries, R and Python packages, etc.) into an [image](https://docs.docker.com/get-started/overview/#images). When you run this image it will give you access to all of the software dependencies of the project inside a completely self-contained instance called a [container](https://docs.docker.com/get-started/overview/#containers). You can then run, explore, and test the project's code inside the container. When you are done, you can stop the container and remove the image from your computer.

To set up the reproducible environment on your computer:

Download this project. Make note of the location you download the project to.

Download and install [Docker Desktop](https://docs.docker.com/get-docker/).

Download the Docker image for this project by running the following command in the terminal:

```
docker pull ghcr.io/mccarthy-m-g/mccarthy-20XX:latest
```

Run the Docker image inside a container:

```
docker run -it -v ~/mccarthy-20XX:/mccarthy-20XX my-docker-image
```

Open [RStudio Server]() in your browser.

Open the RStudio project.

You can now run, explore, and test the project's code.

To remove the reproducible environment from your computer:

Stop the container.

Remove the image.

Delete the project.

## License

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC_BY--NC--ND_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

Except where otherwise noted, this project is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](https://creativecommons.org/licenses/by-nc-nd/4.0/). Please see `LICENSE.md` for a copy of this license.

You are free to:

**Share**---copy and redistribute the material in any medium or format.

Under the following terms:

**Attribution**---You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.

**NonCommercial**---You may not use the material for commercial purposes.

**NoDerivatives**---If you remix, transform, or build upon the material, you may not distribute the modified material.

### Data License

[![License: CC0-1.0](https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg)](http://creativecommons.org/publicdomain/zero/1.0/)

All data inside the `data/EEG-sensor-connectivity/` directory is licensed under the [Creative Commons CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/). Please see `LICENSE.md` for a copy of this license.

You are free to copy, modify, distribute and perform the work, even for commercial purposes, all without asking permission.

### Code License

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

All source code in this project is licensed under the [MIT License](https://opensource.org/licenses/MIT). Please see `LICENSE.md` for a copy of this license.

You are free to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so.

Under the following terms:

**Attribution**---You must include the following copyright notice "Copyright (c) 2022 Michael McCarthy" and provide a link to the license in all copies or substantial portions of the Software.
