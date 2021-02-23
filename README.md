
# TransferTraitDemo
Scripts to use the methodology of **'Linking functional traits to demographic parameters in high-diversity community models'** by L. Chalmandrier, F. Hartig, D.C. Laughlin, H. Lischke, M. Pichler, D.B. Stouffer, L. Pellissier;


The present code allows to calibrate a community model on a mountain grassland dataset via functional traits. The method assumes a linear relationship between functional traits and the demographic parameters of the community model with a transfer function. Through Bayesian inference it then calibrates the transfer function to best reproduce observed species abundances. 

All the necessary code to reproduce the analysis in the main manuscript and in the supporting information (including figures, tables and their associated source data).

## Scripts and functions
### Main analysis
**main/obs.R** is the main script that show the procedure to treat the dataset and reproduces the model used in the original article. Details explanation of each function is provided in their associated script. The results are stored in *results/obs/*

**main/prior.R** contains the specification of the prior. It is directly called by **main/obs.R**.

**main/NM.R** is the script to run the calibration procedure with randomized trait values. The results are stored in **results/NM/Lls.txt**

**main/output_analysis.R** contains a small script to visualize the results of **main/obs.R** that are stored in *results/obs/*

**main/figures.R** contains a script to generate the figures 2 to 4 (stored in ./figures) of the manuscript as well as the source data (stored in the ./SourceData)

**lib/likelihood.R** contains the likelihood function to optimize for the three traits models

**lib/LV_model_wrapped.cpp** contains the ODE model. It is coded in C++ and uses the boost libraries

**lib/trait2demo.R** contains the transfer function that estimates demographic parameters from functional traits.

### Comparative analyses
**main/comparative_analysis/GLM_stan.R** is the script to run the SDM model.

**main/comparative_analysis/LVM_stan.R** is the script to run the latent variable model (jSDM).

**main/comparative_analysis/saturated_stan.R** is the script to run the saturated model.

**main/comparative_analysis/sjSDM_BT.R** is the script to run the sjSDM model with MCMC sampling.

**main/comparative_analysis/sjSDM_MLE.R** is the script to run the sjSDM (jSDM) model. To run the sjSDM with the multinomial response you have to install the multinomial branch:
```r
devtools::install_github("https://github.com/TheoreticalEcology/s-jSDM", ref="multinomial", subdir = "sjSDM")
```
**main/4corner.R** is the script to run the fourth corner analysis and generate Table S6

**main/comparative_analysis/FiguresCA.R** is the script to produce the Supplemenentary figures 6 and 7

**lib/glm.stan** contains the stan model for the SDM.

**lib/lvm.stan** contains the stan model for the latent variable model (jSDM).

**lib/saturated.stan** contains the stan model for the saturated model.

### Supplementary analyses
The folder *main/permutation_analysis/* contains the scripts to generate the analysis of the Supplementary notes 1 that tested if the approach was robust to the shuffling of trait axes. The scripts are written following the structure of the scripts **main/obs.R** and **main/prior.R**. The results of each permutation script are stored in *results/obs_132/*, *results/obs_213/*... depending on the trait order.

The folder *main/four_trait_analysis/* contains the scripts to generate the analysis of the Supplementary notes 1 that tested if the approach was robust to the addition of a fourth axis. The scripts are written following the structure of the scripts **main/obs.R** and **main/prior.R**. The results are stored in *results/obs_n4/* ; the model uses the likelihood function stored in **lib/likelihood_n4.R**

The script **main/figuresSI.R** generates the figures S1-S5 and the table S1

## Data : data/data.Rdata
The data is contained in **data/data.Rdata** and includes:

1. ixp: a matrix with plot as rows and species as columns. The matrice entries contains the number of individuals of each species sampled in each plot.
2. temp: vector of the mean annual temperature in each plot
3. spxt: matrice of the functional traits for each species: reproductive height (Hrepro), vegetative height (Hveg), SLA (specific leaf area), LDMC (leaf dry matter content), LCC (leaf carbon content), LNC (leaf nitrogen content), LD13C (leaf isotopic ratio of carbon 13), LD15N (leaf isotopic ratio of nitrogen 15). 
4. FG : data frame that contains for each species, its full taxonomic name and its associated functional group. Note that three species have ambiguous taxonomic identification and are noted 'carex_sp', 'Ranunculus_sp.", "hieracium sp."

System requirements R 3.6.0
