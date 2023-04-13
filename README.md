# FPCA-SGL

This Github repository contains the R code for the simulation studies in the manuscript titled "Functional principal component analysis and sparse-group LASSO to identify associations between biomarker trajectories and mortality among hospitalized SARS-CoV-2 infected individuals" by Tingyi Cao.

First, the 3 files "Simulation_data_Model1.R", "Simulation_data_Model2.R", and "Simulation_data_Model3.R" are run to simulate data with the 3 models under 5 scenarios (as listed in Table 2 in the manuscript). The simulated datasets are then used in the other 3 files where: 1) the "Simulation_FPCA_SGL.R" file implements the proposed two-stage FPCA-SGL analytic strategy; 2) the "Simulation_baseline_SGL.R" and the "Simulation_peak_SGL.R" files implement the two comparator methods.
