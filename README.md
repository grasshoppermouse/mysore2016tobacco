# mysore2016tobacco

[![DOI](https://zenodo.org/badge/288039354.svg)](https://zenodo.org/badge/latestdoi/288039354)

This repository contains all code necessary to reproduce:

[*The impact of information about tobacco-related reproductive vs. general health risks on South Indian women's tobacco use decisions*](https://grasshoppermouse.github.io/mysore2016tobacco/). Caitlyn D. Placek, Renee E. Magnan, Vijaya Srinivas, Poornima Jaykrishna, Kavitha Ravi, Anisa Khan, Purnima Madhivanan and Edward H. Hagen

Instructions:

1. Clone this repository
2. Open the project in RStudio or `cd` into the directory and launch `R`. This will automatically bootstrap [`renv`](https://rstudio.github.io/renv/index.html).
3. After the bootstrapping process, enter the following in the console: `renv::restore()`. This should install all the necessary packages, including the separate data package [`placektobaccopublic`](https://github.com/grasshoppermouse/placektobaccopublic), in an isolated project-specific library.
4. knit the `paper.Rmd` file using the RStudio GUI or with `rmarkdown::render('paper.Rmd')`. This will generate the preprint file [`paper.html`](https://grasshoppermouse.github.io/mysore2016tobacco/), which will display in the RStudio Viewer or can be viewed in any web browser. (Note: if not using RStudio, you will need a recent version of [pandoc](https://pandoc.org) installed.)

Note 1: Analyses used "R version 3.6.3 (2020-02-29)". You might need to install this version of R to reproduce them.

Note 2: Some results will differ slightly from the official version because identifying information in `placektobaccopublic`, such as age and income, has been aggregated into coarser categories to help preserve anonymity.
