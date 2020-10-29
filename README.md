# tortoise-population-model
A pair of analyses meant to estimate population demographics and project future population viability of Gopher Tortoises in Alabama.

This repository contains five files: 

1) A .csv datafile ("captures.csv") that includes captures of gopher tortoises from various sites in Alabama, including the six study populations in Conecuh National Forest. This is the raw data from which patterns of population demographics are estimated from six sites in Alabama during 1991-2020.

2) A .jags file ("multisite-CMR_fixed-site-effects.jags") that describes a multi-state, multi-site mark-recapture model that can be used to estimate apparent survival, state-transition, and recapture probabilities for tortoises among the six populations. See our manuscript for details about the model and what it does. 

3) A short .R script ("functions.R") that contains two functions useful for fitting the mark-recapture model: i) a function to specify the known latent state of each individual, and ii) a function to create initial values and specify NAs for the first encounter 

4) An R script ("fit-multisite-CMR-fixed-effects.R") that can be used to analyze the capture data (.csv file) with the multistate mark-recapture model (.jags file) to estimate demographics. The "functions.R" file is necessary here also. After running this script, it's important to save the results of the mark-recapture analysis as an .RDS file, to then import into other scripts for visualization with graphs and for downstream analysis. 

5) An R script ("tortoise-PVA-conecuh.R") that can be used to i) visualize the results of the estimation script (output of "fit-multisite-CMR-fixed-effects.R"), ii) perform population projections for each of the study populations using algebraic formulas, and iii) perform elasticity analysis and other analyses for each population. 

I also wrote scripts to simulate data for our multi-state framework; analysis of simulated data suggested that our model is pretty accurate with little bias. If you have questions about these analyses or are interested in the script for the simulated data, feel free to contact me at brian dot folt at gmail dot com. 
