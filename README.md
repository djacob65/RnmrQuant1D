# RnmrQuant1D
Package dedicated to quantification from 1D NMR spectra, including peak fitting and calibration using standard spectra.

Visit the [Wiki](https://github.com/djacob65/RnmrQuant1D/wiki/home/) page for a detailed tutorial

------

### Installation

1. **R version >= 4.3.0** - [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools44/rtools.html) >=4.3 required for Windows - Recommanded : [RStudio](https://posit.co/download/rstudio-desktop/)

2. **Rnmr1D (>= 1.6.0)**

* Packages to be installed : (see https://github.com/INRA/Rnmr1D)
  ```r
  Rcpp (>= 0.12.7), base64enc (>= 0.1), MASS(>= 7.3), Matrix, methods, scales, doParallel (>= 1.0.11), 
  foreach (>= 1.4.4), igraph (>= 1.2.1), impute (>= 1.54.0)(*) , MassSpecWavelet (>= 1.46.0)(*), 
  ptw (>= 1.9), signal (>= 0.7), XML (>= 3.98), ggplot2 (>= 3.0.0), plotly (>= 4.8.0), plyr (>= 1.8.4), 
  minqa (>= 1.2.4)
  ```
    * (*) packages deposited in Bioconductor

3. **RnmrQuant1D**

* Packages to be installed
  ```r
  Rnmr1D (>= 1.6.0), foreach (>= 1.5.2), doParallel (>= 1.0.17), magrittr (>= 2.0.3), 
  DT (>= 0.29), plotly (>= 4.8.0), openxlsx (>= 4.2.0),  webshot2 (>= 0.1.1)
  ```

* Installation of the R package
  ```r
  if (!require("devtools"))
    install.packages("devtools", repos="https://cran.rstudio.com")
  devtools::install_github("djacob65/RnmrQuant1D", dependencies = TRUE)

  ```

------

### Funded by:

* Agence Nationale de la Recherche - [ANR-21-CE21-0014](https://anr.fr/Project-ANR-21-CE21-0014)
* [UMR INRAE BFP Bordeaux](https://eng-bfp.bordeaux-aquitaine.hub.inrae.fr/)
* [UR INRAE BIA Nantes](https://eng-ur-bia.angers-nantes.hub.inrae.fr/)

### License

GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007 - See the included LICENSE file.
