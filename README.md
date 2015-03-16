# Replication material for: Deterministic and stochastic trends in the Lee-Carter mortality model
### Laurent Callot, Niels Haldrup, Malene Kallestrup Lamb

---

Author: Laurent Callot (l.callot@vu.nl)

Date: 19/11/2014


This repository contains the replication material for _Deterministic and stochastic trends in the Lee-Carter mortality model_. All the computations are carried using *R* and the *knitr* pacakge for easy replication. The material is divided between 4 folders:

 - results contains the .Rnw file used to generate the output in the paper. To replicate the results run R, load the knitr package and simply type: _knit2pdf(pubplots.Rnw)_ to generate a pdf with the plots and tables of the paper.  
 - code contains some of the functions used to carry the computations.
 - logm contains the smoothed log mortality rates in csv format. It is divided in 4 folders containing the data at different starting dates.
 - hdm_data contains the code used to transform the raw data from the human mortality database to log mortality rates. The original data cannot be provided due to the license agreement of HMD. 
  
