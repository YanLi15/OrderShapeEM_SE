# OrderShapeEM_SE: Identifying SE genes from spatially resolved transcriptomic data using OrderShapeEM

## Overview

This repository contains the R implementation of simulations and case studies as well as underlying data in our manuscript: Covariate-assistant statistical analysis of spatially resolved transcriptomic data incorporating multi-omics auxiliary information. 

All the functions for data analysis are available in ‘R’ folder. The results of real data guided simulations in different settings and case studies are saved in ‘Simulation’ and ‘Output’ folders. The R codes and scripts used to reproduce the results can be found in ‘Analysis ’ folder. The datasets downloaded from publicly available databases for validation of our experiments are also included in ‘ValidateData’ folder. The raw count spatial transcriptomic data and auxiliary data used in the case studies can be downloaded from the links provided in the “Data and codes availability ” Section of the manuscript and stored in ‘RealData’ and ‘Covariates’ folders, respectively. 

The original [SPARK](https://github.com/xzhoulab/SPARK ) and [OrderShapeEM](https://github.com/jchen1981/OrderShapeEM) are required for the integrative analysis.

## Work flow

We put  three case studies and a newly added case study in analyzing spatially resolved transcriptomic data into the ‘Analysis’ folder, including the Mouse olfactory bulb (MOB) dataset, Mouse cerebellum (MC) dataset, human reast cancer (HBC) dataset and human heart (HH) dataset. Taking the MOB data as an example, the integrative analysis includes the following steps:

- MOB\_data.R: load the spatial transcriptomic data and analyze the data with SPARK; load the auxiliary data, i.e., single-cell RNA-seq (scRNA-seq) data, from MOB and make differential analysis between different cell-types with Seurat.  The results of SPARK and scRNA-seq analysis are saved as ‘MOB.RData’ in the ‘Output’ folder.
- MOB\_analysis.R: load the data stored in ‘MOB.RData’, take $p$-values obtained from SPARK as primary $p$-values, $p$-values obtained from scRNA-seq data as auxiliary $p$-values. Use Cauchy combination rule to combine $p$ values from multiple studies into one covariate. After matching them by gene, run BY, BH, SABHA and OrderShapeEM procedures to identify statistically significant spatially expressed genes with multiple testing adjustment. The SE genes identified by different methods are returned and saved as ‘MOB_result.RData’ in the ‘Output’ folder.
- MOB\_plot.R: summarize and visualize the spatial expression patterns using the hierarchical clustering analysis based on the SE genes identified by OrderShapeEM.
- MOB_plot.R: visualize the spatial expression patterns of specific genes.
- MOB\_validate.R: validate the SE genes identified by OrderShapeEM with three additional lines of evidence from other published literature or databases.
- MOB\_go.R: perform GO enrichment analysis on SE genes detected by OrderShapeEM.

To ease comparison, we code all other comparing methods in addition to OrderShapeEM in the RunAll() function. We construct a Cauchy() function to code the Cauchy combination rule for auxiliary data from multi-omics studies. An R script named ‘example.R’ is added to showcase the work flow of our integrative analysis method (loading data, SPARK, Cauchy combination, OrderShapeEM).

*Note: For large-scale dataset, we have made a small change in the spark.test() function from R package SPARK to close the connections in time to avoid overload errors. In the analysis of mouse cerebellum and human heart data using SPARK, we re-sourced the modified spark.test() function stored in the ‘R’ folder to make the computation feasible.*

