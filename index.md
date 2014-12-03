---
layout: page
title: quantro
tagline: Additional Materials and Scripts
---

The additional materials (R code and scripts) for the paper *When to use Quantile Normalization?* authored by Stephanie C. Hicks and Rafael A. Irizarry are provided here. 

## R-packages

#### Supplementary Sections 1 and 3

These are the static R-packages (`quantro` and `quantroSim`) used in the analysis for the paper (Supplementary Material Sections 1 and 3, respectively). 

* [quantro_0.99.1.tar.gz](https://github.com/stephaniehicks/quantroPaper/raw/master/Rpkgs/quantro_0.99.1.tar.gz): Most recent version available on [Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/quantro.html)
* [quantroSim_0.0.1.tar.gz](https://github.com/stephaniehicks/quantroPaper/raw/master/Rpkgs/quantroSim_0.0.1.tar.gz): Most recent version available on [Github](https://github.com/stephaniehicks/quantroSim)


## Scripts

#### Supplementary Section 2

In this section, we apply `quantro` to several gene expression and DNA methylation datasets. For a description of all the data used, see Table 1 in the [Supplementary Material](). 

Gene expression (RNA-Seq)

* [pickrellRNASeq](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/geneExpression/pickrellRNASeq.Rmd)
* [mouseStrainsRNASeq](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/geneExpression/mouseStrainsRNASeq.Rmd)

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
* [mycAffyData](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/geneExpression/mycAffyData.Rmd)

DNA methylation (Microarrays)

* [adiposeExerciseMethyl](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/dnaMethylation/adiposeExerciseMethyl.Rmd)
* [pancreaticT2DMethyl]()
* [cellcompMethyl]()



#### Supplementary Section 4

In this section, we perform a simulation study to assess the performance of our method in the `quantro` R/Bioconductor package using the `quantroSim` R-package. 

* [Script for Supplementary Figures 21-22](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/pDiffRandom.R): Bias and MSE
* [Script for Supplementary Figure 23](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/FDR.R): False discovery plots
* [Script for Supplementary Figure 24](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/ROC.R): ROC curves
* [Helper functions](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/quantro-functions.R)