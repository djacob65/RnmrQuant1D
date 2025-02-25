# RnmrQuant1D
Package dedicated to quantification from 1D NMR spectra, including peak fitting and calibration using standard spectra.

<br>

![image](https://github.com/user-attachments/assets/8b19fe92-82a6-4105-bc99-19d6122f16e5)

<br>

See the [Wiki](https://github.com/djacob65/RnmrQuant1D/wiki/home/) page for more details


### Installation

1. **R version >= 4.3.0** - Rtools >=4.3 required for Windows - Recommanded : RStudio

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
  Rnmr1D (>= 1.6.0), openxlsx (>= 4.2.0), data.table (>= 1.15.2), reshape2 (>= 1.4.4),  digest (>= 0.6.35), 
  jsonlite (>= 1.8.8),  repr (>= 1.1.6), webshot2 (>= 0.1.1),  DT (>= 0.29), repr (>= 1.1.6), magrittr (>= 2.0.3)
  ```

* Installation of the R package
  ```r
  require(devtools)
  install_github("djacob65/RnmrQuant1D", dependencies = TRUE)

  ```

### Note

The source code is currently hosted on another repository with private access for sharing to authorized members. When the embargo period is over, the source code will be available here. Thank you for your understanding.

------

### Funded by:

* Agence Nationale de la Recherche - [ANR-21-CE21-0014](https://anr.fr/Project-ANR-21-CE21-0014)
* UMR INRAE BFP Bordeaux
* UR INRAE BIA Nantes

### License

GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007 - See the included LICENSE file.

