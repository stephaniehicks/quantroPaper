---
layout: page
title: quantro
tagline: Additional Materials and Scripts
---

The additional materials (R code and scripts) for the paper *When to use Quantile Normalization?* authored by Stephanie C. Hicks and Rafael A. Irizarry are provided here. 

# R-packages

These are the R-packages (static) used in the analysis for the paper (Supplementary Material Sections 1 and 3). 

* [quantro_0.99.1.tar.gz](): Most recent version available on [Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/quantro.html)
* [quantroSim_0.0.1.tar.gz](): Most recent version available on [Github](https://github.com/stephaniehicks/quantroSim)


# Scripts

#### Supplementary Section 2: Applying `quantro` to gene expression and DNA methylation data

For a description of all the data used, see Table 1 in the [Supplementary Material](). 

Gene expression (RNA-Seq)

* [pickrellRNASeq]()
* [mouseStrainsRNASeq]()

Gene expression (Microarrays)

* [alveolarSmokingAffyData]()
* [lungCOPDAffyData]()
* [brainParkinsonsAffyData]()
* [brainLiverAffyData]()
* [lungCancerAffyData]()
* [breastCancerAffyData]()
* [prostateCancerAffyData]()
* [thyroidCancerAffyData]()
* [stomachCancerAffyData]()
* [liverCancerAffyData]()
* [liverNAFLDAffyData]()
* [mycAffyData]()

DNA methylation (Microarrays)

* [adiposeExerciseMethyl]()
* [pancreaticT2DMethyl]()
* [cellcompMethyl]()



#### Supplementary Section 4: Simulation Study

* [Script for Supplementary Figures 21-22](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/pDiffRandom.R): Bias and MSE
* [Script for Supplementary Figure 23](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/FDR.R): False discovery plots
* [Script for Supplementary Figure 24](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/ROC.R): ROC curves
* [Helper functions](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/quantro-functions.R)