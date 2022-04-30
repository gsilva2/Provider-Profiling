library(MASS)
library(boot)

############## Introduction ##############

# The scripts included in this GitHub are not identical to the analyses included in the manuscript
# In the manuscript, actual patient-level provider data was used for both the simulation and the observed outcome analysis
# In this illustration, we simulate patient-level provider data to illustrate how the methodology and code in the manuscript can be implemented
# After running Script 1, Section A, the resulting data files ("data_prov_simulated.csv", "data_cont_simulated.csv", "data_cat_simulated.csv") can be used to cluster the J (or n_prov) providers into M groups 
# The clustering occurs when script2_simdata_fitn_2cont2cat.py is run. This script has been prepared for M=3 but can be updated for other values of M
# After running the clustering script across different values of M and obtaining the posterior draws for "category" for each value of M, these can be used to select the optimal value of M, M*.
# Selecting the optimal value of M* using the average silhouette method is illustrated in script3_clusteranalysis.R
# Once M* has been identified, the outcome data simulated in Script 1, Section B (using a similar set-up as to what was described in the manuscript), can be used to understand the impact of clustering 
# script4_outcomeanalysis.R shows how the estimation of the SRR within each of the M* groups compares to the estimation of the SRR with a single group and how these can impact the false positive rate (FPR), the false negative rate (FNR), and accuracy

# Order of Analyses:
# First run script1_simulatesampledata.R to generate simulated covariate data (Section A) and outcome data (Section B)
# Then run script2_simdata_fitn_2cont2cat.py for different values of M to cluster the providers using the generated covariate data 
# Then run script3_clusteranalysis.R to process the cluster results for different values of M and select M*
# Then run script4_outcomeanalysis.R to complete analyses with outcome data (Section C) and understand the impact of clustering

############## Section A: Simulating patient covariate data for clustering of providers ##############

# Suppose J=number of providers=n_prov is 120, similar to application
n_prov<-120
  
# Generating number of patients per provider
set.seed(12345)
patients_per_prov<-round(c(rnorm(n_prov-30, 80, sd=20), rnorm(30, 320, sd=50)))
total_patients<-sum(patients_per_prov)
  
# Suppose data is being generated from 3 different clusters
cluster_prov<-c(rep(1, n_prov/3), rep(2, n_prov/3), rep(3, n_prov/3))

Provider<-rep(1:120, times=patients_per_prov)
Cluster<-rep(cluster_prov, times=patients_per_prov)

mu_1<-c(30, 6)
mu_2<-c(40, 3)
mu_3<-c(45, 8)

X_cont<-matrix(NA, ncol=2, nrow=total_patients)
for(i in 1:total_patients){
  set.seed(i*10)
  if(Cluster[i]==1){X_cont[i, ]<-mvrnorm(n = 1, mu_1, Sigma=matrix(c(5, 0.25, 0.25, 1.5), nrow=2), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)}
  if(Cluster[i]==2){X_cont[i, ]<-mvrnorm(n = 1, mu_2, Sigma=matrix(c(5, 0.25, 0.25, 1.5), nrow=2), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)}
  if(Cluster[i]==3){X_cont[i, ]<-mvrnorm(n = 1, mu_3, Sigma=matrix(c(5, 0.25, 0.25, 1.5), nrow=2), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)}
}

# Scaling continuous covariates for analysis
X_cont[,1]<-scale(X_cont[,1])
X_cont[,2]<-scale(X_cont[,2])

X_cat<-matrix(NA, ncol=4, nrow=total_patients)
for(i in 1:total_patients){
  set.seed(i*10)
  if(Cluster[i]==1){X_cat[i, ]<-as.vector(rmultinom(1, 1, c(0.1, 0.2, 0.6, 0.1)))}
  if(Cluster[i]==2){X_cat[i, ]<-as.vector(rmultinom(1, 1, c(0.6, 0.2, 0.1, 0.1)))}
  if(Cluster[i]==3){X_cat[i, ]<-as.vector(rmultinom(1, 1, c(0.1, 0.4, 0.1, 0.4)))}
}

# Save data files to run script2_simdata_fitn_2cont2cat.py
#write.table(Provider-1, "data_prov_simulated.csv", row.names=FALSE, col.names=FALSE, sep=",")
#write.table(X_cont, "data_cont_simulated.csv", row.names=FALSE, col.names=FALSE, sep=",")
#write.table(X_cat, "data_cat_simulated.csv", row.names=FALSE, col.names=FALSE, sep=",")

############## Section B: Simulating patient outcome data for analysis of outcome once providers are clustered ##############

# Changing format of categorical variables for further analyses 
categorical_1<-X_cat[,1]+X_cat[,2]
categorical_2<-X_cat[,1]+X_cat[,3]

# Creating a single dataset with all the information
data<-cbind(Provider, X_cont, categorical_1, categorical_2)
colnames(data)<-c("Provider", "age", "hosp_los", "female", "renal")
data<-as.data.frame(data)
data$Provider<-as.factor(data$Provider)

######### Generating Probabilities for Binary Outcomes ######### 

###### Parameter specification ###### 
# These match what was specified in Section 4 of paper

# Parameters for linear outcome model
mean_j<--1.75
beta_1<-0.020
beta_2<-0.30
beta_3<--0.25
beta_4<-0.45
beta_5<-0.0


# For this simulation we will assume there are outliers
# To simulate no outliers, tau_j2 would be replaced with 0
tau_j2<-(0.1)^2

set.seed(5)
alpha_j1<-rnorm(n_prov, mean_j, sd=sqrt(tau_j2))
alpha_j2<-rnorm(n_prov, 1+mean_j, sd=sqrt(tau_j2))

# In the manuscript, the percent of outliers examined was 10, 15, 20% 
# Suppose the percent of outliers is 15% 
# This means there will be 6 providers in each cluster who are outliers 0.15*n_prov=18
# The last 6 providers in each cluster are outliers
w<-rep(c(rep(0, 34), rep(1, 6)), length=120)
# Saving which providers are outliers
outliers<-w>0
#save(outliers, file="outliers_sim_15.RData")

# Providers who are outliers will have alpha_j2, providers who are not outliers will have alpha_j1 
alpha_j<-w*alpha_j2+(1-w)*alpha_j1


# Now that each of the n_prov providers has an alpha_j attached to them, remap each providers alpha_j to their corresponding patients
hi<-t(rbind(1:n_prov, alpha_j))
colnames(hi)<-c("Provider", "intercept")
data<-merge(data, hi, by="Provider")
alpha_j<-data$intercept


###### Probabilities of the outcome under different links, model specifications ###### 

# Correct link, correct model 
model<-alpha_j+beta_1*data$age+beta_2*data$hosp_los+beta_3*data$female+beta_4*data$renal
log_prop<-inv.logit(model)
# Incorrect link, correct model
box1_prop<-1-pnorm(-(model))^0.5

beta_5<-1.0
beta_6<-1.0
model<-alpha_j+beta_1*data$age+beta_2*data$hosp_los+beta_5*(data$hosp_los)^2+beta_6*(data$age)*(data$hosp_los)+beta_3*data$female+beta_4*data$renal
# Correct link, incorrect model 
log_imp1<-inv.logit(model)
# Incorrect link, incorrect model
box1_imp1<-1-pnorm(-(model))^0.5


beta_5<-1.5
beta_6<-1.5
model<-alpha_j+beta_1*data$age+beta_2*data$hosp_los+beta_5*(data$hosp_los)^2+beta_6*(data$age)*(data$hosp_los)+beta_3*data$female+beta_4*data$renal
# Correct link, incorrect model 
log_imp2<-inv.logit(model)
# Incorrect link, incorrect model
box1_imp2<-1-pnorm(-(model))^0.5



###### Function to obtain simulated outcomes and full dataset for future analyses ###### 
# Number of outcome draws
rep=200

data_generator<-function(probab){
  data_compile<-list()
  for(j in 1:rep){
    set.seed(10*j)
    data$sim_Y<-rbinom(dim(data)[1], 1, probab)
    data_compile[[j]]<-data
  }
  return(data_compile)
}

# Here we illustrate with log_prop (proper link, proper model), log_imp1 (proper link, incorrect model)

log_prop<-data_generator(log_prop)
#box1_prop<-data_generator(box1_prop)

log_imp1<-data_generator(log_imp1)
#box1_imp1<-data_generator(box1_imp1)

#log_imp2<-data_generator(log_imp2)
#box1_imp2<-data_generator(box1_imp2)

# Save list of datasets for SRR, FPR/FNR calculations
#save(log_prop, file="logprop_dat.Rdata")
#save(box1_prop, file="box1prop_dat.Rdata")
#save(log_imp1, file="logimp1_dat.Rdata")
#save(box1_imp1, file="box1imp1_dat.Rdata")
#save(log_imp2, file="logimp2_dat.Rdata")
#save(box1_imp2, file="box1imp2_dat.Rdata")

