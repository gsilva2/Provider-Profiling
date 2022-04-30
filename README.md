# Provider-Profiling
# Introduction 

### The scripts included in this GitHub are not identical to the analyses included in the manuscript
### In the manuscript, actual patient-level provider data was used for both the simulation and the observed outcome analysis
### In this illustration, we simulate patient-level provider data to illustrate how the methodology and code in the manuscript can be implemented
### After running Script 1, Section A, the resulting data files ("data_prov_simulated.csv", "data_cont_simulated.csv", "data_cat_simulated.csv") can be used to cluster the J (or n_prov) providers into M groups 
### The clustering occurs when script2_simdata_fitn_2cont2cat.py is run. This script has been prepared for M=3 but can be updated for other values of M
### After running the clustering script across different values of M and obtaining the posterior draws for "category" for each value of M, these can be used to select the optimal value of M, M*.
### Selecting the optimal value of M* using the average silhouette method is illustrated in script3_clusteranalysis.R
### Once M* has been identified, the outcome data simulated in Script 1, Section B (using a similar set-up as to what was described in the manuscript), can be used to understand the impact of clustering 
### script4_outcomeanalysis.R shows how the estimation of the SRR within each of the M* groups compares to the estimation of the SRR with a single group and how these can impact the false positive rate (FPR), the false negative rate (FNR), and accuracy

### Order of Analyses:
### First run script1_simulatesampledata.R to generate simulated covariate data (Section A) and outcome data (Section B)
### Then run script2_simdata_fitn_2cont2cat.py for different values of M to cluster the providers using the generated covariate data 
### Then run script3_clusteranalysis.R to process the cluster results for different values of M and select M*
### Then run script4_outcomeanalysis.R to complete analyses with outcome data (Section C) and understand the impact of clustering
