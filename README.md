These materials are a computationally reproducible version of the paper:

Smith, C. & Baker, D.H. (2025). Neural correlates of the deployment of spatial attention, and their modulation by repetitive movements. *PLoS ONE*, 20(9): e0332583. https://doi.org/10.1371/journal.pone.0332583

The file manuscript.qmd is a Quarto markdown file that will perform all analyses and figure creation, and produce a pdf version of the manuscript.

The full repository can be downloaded (cloned), and will attempt to download any data files required from the OSF repository for this project:
http://doi.org/10.17605/OSF.IO/YUHZ4

The 'docker' directory contains a Dockerfile and instructions for making a local computationally reproducible version of the analysis. In addition, the Docker environment is set up to run automatically on a remote server via Github Actions, each time a change is made (i.e. on a 'commit' to the repo). The output document is then posted back to the main repository (manuscript.pdf). If you want to make changes to the analysis and have these build automatically, you can fork the repository into your own account.

![autobuild](https://github.com/bakerdh/attentionstimming/workflows/autobuild/badge.svg)