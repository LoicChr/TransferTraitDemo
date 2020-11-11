# TransferTraitDemo
Scripts to use the methodology of 'Linking functional traits to demographic parameters in high-diversity community models' by Loïc Chalmandrier*, Florian Hartig, Daniel Laughlin, Heike Lischke, Daniel Stouffer & Loïc Pellissier
* main code contributor

The present code allows to calibrate a community model on a mountain grassland dataset via functional traits. The method assumes a linear relationship between functional traits and the demographic parameters of the community model with a transfer function. Through Bayesian inference it then calibrates the transfer function to best reproduce observed species abundances.  

## Scripts and functions
main/obs.R is the main script that show the procedure to treat the dataset and reproduces the model used in the original article. Details explanation of each function is provided in their associated script.
main/prior.R contains the specification of the prior. It is directly called by main/obs.R
main/output_analysis.R contains a small script to visualize the results of main/obs.R that are stored in results/obs/

lib/likelihood.R contains the likelihood function to optimize
lib/LV_model_wrapped.cpp contains the ODE model. It is coded in C++ and uses the boost libraries
lib/trait2demo.R contains the transfer function that estimates demographic parameters from functional traits.

## Data file : data/data.Rdata
This Rdata file contains the following objects
- ixp: a matrix with plot as rows and species as columns. The matrice entries contains the number of individuals of each species sampled in each plot.
- temp: vector of the mean annual temperature in each plot
- spxt: matrice of the functional traits for each species: reproductive height (Hrepro), vegetative height (Hveg), SLA (specific leaf area), LDMC (leaf dry matter content), LCC (leaf carbon content), LNC (leaf nitrogen content), LD13C (leaf isotopic ratio of carbon 13), LD15N (leaf isotopic ratio of nitrogen 15). 
- FG : data frame that contains for each species, its full taxonomic name and its associated functional group. Note that three species have ambiguous taxonomic identification and are noted 'carex_sp', 'Ranunculus_sp.", "hieracium sp."

## Result file: results/obs/*
This folder contains eight chains stored in individual Rdata files. Those are the chains run in parallel with main/obs.R that were used to in the original manuscript. They are visualized using main/output_analysis.R

System requirements R 3.6.0
