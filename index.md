---
layout: page
title: quantro
tagline: Additional Materials and Scripts
---

The additional materials (R code and scripts) for the paper *When to use Quantile Normalization?* authored by Stephanie C. Hicks and Rafael A. Irizarry are provided here. 

## R-packages

#### quantro

The `quantro` R/Bioconductor package is discussed in detail in the Supplementary Section 1 of the paper.   This package contains a data-driven test for global differences between groups of distributions which asses whether global normalization methods such as quantile normalization should be applied. 

* [quantro_0.99.1.tar.gz](https://github.com/stephaniehicks/quantroPaper/raw/master/Rpkgs/quantro_0.99.1.tar.gz) (static)
* [quantro vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/quantro/inst/doc/quantro-vignette.pdf)
* Most recent version available on [Bioconductor](http://www.bioconductor.org/packages/release/bioc/html/quantro.html)


#### quantroSim 

The `quantroSim` R package is discussed in detail in the Supplementary Section 3 of the paper.  This is the supporting data simulation package for the `quantro` R/Bioconductor package which can be used to simulate gene expression and DNA methylation data. 

* [quantroSim_0.0.1.tar.gz](https://github.com/stephaniehicks/quantroPaper/raw/master/Rpkgs/quantroSim_0.0.1.tar.gz) (static)
* [quantroSim vignette](https://github.com/stephaniehicks/quantroSim/raw/master/vignettes/quantroSim-vignette.pdf)
* Most recent version available on [Github](https://github.com/stephaniehicks/quantroSim)



## Scripts

#### Applications to gene expression and DNA methylation

In this section, we apply `quantro` to several gene expression and DNA methylation datasets. For a description of all the data used, see Table 1 in the [Supplementary Material]().  Below are Rmarkdown files containing the analysis using `quantro`. 

Gene expression (RNA-Seq)

* [pickrellRNASeq](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/geneExpression/pickrellRNASeq.Rmd)
* [mouseStrainsRNASeq](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/geneExpression/mouseStrainsRNASeq.Rmd)

Gene expression (Microarrays)

* [alveolarSmokingAffyData](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/geneExpression/alveolarSmokingAffyData.Rmd)
* [lungCOPDAffyData](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/geneExpression/lungCOPDAffyData.Rmd)
* [brainParkinsonsAffyData](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/geneExpression/brainParkinsonsAffyData.Rmd)
* [brainLiverAffyData]()
* [lungCancerAffyData]()
* [breastCancerAffyData]()
* [prostateCancerAffyData]()
* [thyroidCancerAffyData]()
* [stomachCancerAffyData]()
* [liverCancerAffyData]()
* [liverNAFLDAffyData](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/geneExpression/liverNAFLDAffyData.Rmd)
* [mycAffyData](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/geneExpression/mycAffyData.Rmd)

DNA methylation (Microarrays)

* [adiposeExerciseMethyl](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/dnaMethylation/adiposeExerciseMethyl.Rmd)
* [pancreaticT2DMethyl](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/dnaMethylation/pancreaticT2DMethyl.Rmd)
* [cellcompMethyl](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/dnaMethylation/cellcompMethyl.Rmd)



#### Simulation studies

In this section, we perform several Monte Carlo simulation studies to assess the performance of our method in the `quantro` R/Bioconductor package using the `quantroSim` R-package (discussed in further detail in Supplementary Section 4). 

* [Script for Supplementary Figures 21-22](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/pDiffRandom.R): Bias and MSE
* [Script for Supplementary Figure 23](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/FDR.R): False discovery plots
* [Script for Supplementary Figure 24](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/ROC.R): ROC curves
* [Helper functions needed for simulation studies](https://github.com/stephaniehicks/quantroPaper/blob/master/scripts/quantroSimStudy/quantro-functions.R)