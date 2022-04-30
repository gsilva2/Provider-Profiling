import numpy as np
import pymc as mc
import random
import math as mt
from numpy import genfromtxt
from functools import reduce
from operator import add

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


#################################### Importing Data ####################################

# Make sure continuous data file and categorical data file have no column or row names
# Categorical data file should have number of rows equal to the number of patients/units, number of columns equal to the product of the number of levels across all categorical covariates, D

data = genfromtxt('data_cont_simulated.csv', delimiter=',')
data=data.astype(np.float64)
data_cat= genfromtxt('data_cat_simulated.csv', delimiter=',')
data_cat=np.array(data_cat, dtype='int')
# This vector will specify which cell of the contingency table with D potential entries that element belongs to
position=(data_cat!=0).argmax(axis=1)


prov1= genfromtxt('data_prov_simulated.csv', delimiter=',')
prov1=np.array(prov1, dtype='int')


#################################### Model Specification and Fitting ####################################

# Truth in simulated dataset is three clusters
# To fit other number of clusters, simply change the value of n to the required value and then uncomment the parameters needed in the model specification section

#n=2
n=3
#n=4
#n=5
#n=6

# Number of continuous covariates
p = data.shape[1] 	
# Number of potential categories among categorical covariates (equal to the product of the number of levels across all categorical covariates)
D = data_cat.shape[1]  
nprov=len(np.unique(prov1))


################# Specification of Parameters, Priors ################# 

####### Mean #######

prior_mu=np.zeros(p)

mu1_1 = mc.Normal("means1_1",prior_mu,.01,size=p)
mu1_2 = mc.Normal("means1_2",prior_mu,.01,size=p)
mu1_3 = mc.Normal("means1_3",prior_mu,.01,size=p)
mu1_4 = mc.Normal("means1_4",prior_mu,.01,size=p)
mu_1=np.array([mu1_1, mu1_2, mu1_3, mu1_4], object)

mu2_1 = mc.Normal("means2_1",prior_mu,.01,size=p)
mu2_2 = mc.Normal("means2_2",prior_mu,.01,size=p)
mu2_3 = mc.Normal("means2_3",prior_mu,.01,size=p)
mu2_4 = mc.Normal("means2_4",prior_mu,.01,size=p)
mu_2=np.array([mu2_1, mu2_2, mu2_3, mu2_4], object)

mu3_1 = mc.Normal("means3_1",prior_mu,.01,size=p)
mu3_2 = mc.Normal("means3_2",prior_mu,.01,size=p)
mu3_3 = mc.Normal("means3_3",prior_mu,.01,size=p)
mu3_4 = mc.Normal("means3_4",prior_mu,.01,size=p)
mu_3=np.array([mu3_1, mu3_2, mu3_3, mu3_4], object)

#mu4_1 = mc.Normal("means4_1",prior_mu,.01,size=p)
#mu4_2 = mc.Normal("means4_2",prior_mu,.01,size=p)
#mu4_3 = mc.Normal("means4_3",prior_mu,.01,size=p)
#mu4_4 = mc.Normal("means4_4",prior_mu,.01,size=p)
#mu_4=np.array([mu4_1, mu4_2, mu4_3, mu4_4], object)

#mu5_1 = mc.Normal("means5_1",prior_mu,.01,size=p)
#mu5_2 = mc.Normal("means5_2",prior_mu,.01,size=p)
#mu5_3 = mc.Normal("means5_3",prior_mu,.01,size=p)
#mu5_4 = mc.Normal("means5_4",prior_mu,.01,size=p)
#mu_5=np.array([mu5_1, mu5_2, mu5_3, mu5_4], object)

#mu6_1 = mc.Normal("means6_1",prior_mu,.01,size=p)
#mu6_2 = mc.Normal("means6_2",prior_mu,.01,size=p)
#mu6_3 = mc.Normal("means6_3",prior_mu,.01,size=p)
#mu6_4 = mc.Normal("means6_4",prior_mu,.01,size=p)
#mu_6=np.array([mu6_1, mu6_2, mu6_3, mu6_4], object)

#meansAll= np.array([mu_1, mu_2], object)
meansAll= np.array([mu_1, mu_2, mu_3], object)
#meansAll= np.array([mu_1, mu_2, mu_3, mu_4], object)
#meansAll= np.array([mu_1, mu_2, mu_3, mu_4, mu_5], object)
#meansAll= np.array([mu_1, mu_2, mu_3, mu_4, mu_5, mu_6], object)


####### Precision #######

prec_mat_1 = mc.Wishart("prec_matrix_1",n,np.eye(p)*.1)
prec_mat_2= mc.Wishart("prec_matrix_2",n,np.eye(p)*.1)
prec_mat_3 = mc.Wishart("prec_matrix_3",n,np.eye(p)*.1)
#prec_mat_4 = mc.Wishart("prec_matrix_4",n,np.eye(p)*.1)
#prec_mat_5 = mc.Wishart("prec_matrix_5",n,np.eye(p)*.1)
#prec_mat_6= mc.Wishart("prec_matrix_6",n,np.eye(p)*.1)


#precsAll= np.array([prec_mat_1, prec_mat_2], object)
precsAll= np.array([prec_mat_1, prec_mat_2, prec_mat_3], object)
#precsAll= np.array([prec_mat_1, prec_mat_2, prec_mat_3, prec_mat_4], object)
#precsAll= np.array([prec_mat_1, prec_mat_2, prec_mat_3, prec_mat_4, prec_mat_5], object)
#precsAll= np.array([prec_mat_1, prec_mat_2, prec_mat_3, prec_mat_4, prec_mat_5, prec_mat_6], object)


####### Category Probabilities #######

pi_1=mc.Dirichlet('pi_m1', theta=(1,)*D)
pi_2=mc.Dirichlet('pi_m2', theta=(1,)*D)
pi_3=mc.Dirichlet('pi_m3', theta=(1,)*D)
#pi_4=mc.Dirichlet('pi_m4', theta=(1,)*D)
#pi_5=mc.Dirichlet('pi_m5', theta=(1,)*D)
#pi_6=mc.Dirichlet('pi_m6', theta=(1,)*D)


#piAll=np.array([pi_1, pi_2], object)
piAll=np.array([pi_1, pi_2, pi_3], object)
#piAll=np.array([pi_1, pi_2, pi_3, pi_4], object)
#piAll=np.array([pi_1, pi_2, pi_3, pi_4, pi_5], object)
#piAll=np.array([pi_1, pi_2, pi_3, pi_4, pi_5, pi_6], object)


####### Latent Probabilities #######

pip = mc.Dirichlet('pip', theta=(1,)*n)
category = mc.Categorical('category', p=pip, size=nprov)


################# Specification of Deterministic Variables #################

@mc.deterministic
def mean(category=category, meansAll=meansAll):
   lat = category[prov1]
   new = meansAll[lat, position]
   return new
 
@mc.deterministic
def prec(category=category, precsAll=precsAll):
       lat = category[prov1]
       return precsAll[lat]


@mc.deterministic
def pi1(category=category, piAll=piAll):
    lat=category[prov1]
    ans=piAll[lat]
    return ans
   

################# Observed Likelihood #################
    
@mc.observed
def obs(value=data, mean=mean, prec=prec, pi1=pi1,  position=position):
   return sum(mc.mv_normal_like(v, m, T)+float(np.log(np.append(prob1, 1-sum(prob1))[pos])) for v,m,T,prob1,pos in zip(data, mean, prec, pi1, position))


################# Specifying Starting Values #################


start=np.tile(range(0, n), mt.floor(nprov/n))
category.value=np.array(start)

mu1_1.value=[0, 0]
mu1_2.value=[0, 0]
mu1_3.value=[0, 0]
mu1_4.value=[0, 0]
mu2_1.value=[0, 0]
mu2_2.value=[0, 0]
mu2_3.value=[0, 0]
mu2_4.value=[0, 0]
mu3_1.value=[0, 0]
mu3_2.value=[0, 0]
mu3_3.value=[0, 0]
mu3_4.value=[0, 0]
#mu4_1.value=[0, 0]
#mu4_2.value=[0, 0]
#mu4_3.value=[0, 0]
#mu4_4.value=[0, 0]
#mu5_1.value=[0, 0]
#mu5_2.value=[0, 0]
#mu5_3.value=[0, 0]
#mu5_4.value=[0, 0]
#mu6_1.value=[0, 0]
#mu6_2.value=[0, 0]
#mu6_3.value=[0, 0]
#mu6_4.value=[0, 0]


cov_matrix_inv1.value=np.eye(p)
cov_matrix_inv2.value=np.eye(p)
cov_matrix_inv3.value=np.eye(p)
#cov_matrix_inv4.value=np.eye(p)
#cov_matrix_inv5.value=np.eye(p)
#cov_matrix_inv6.value=np.eye(p)


pi_1.value=[0.25, 0.25, 0.25]
pi_2.value=[0.25, 0.25, 0.25]
pi_3.value=[0.25, 0.25, 0.25]
#pi_4.value=[0.25, 0.25, 0.25]
#pi_5.value=[0.25, 0.25, 0.25]
#pi_6.value=[0.25, 0.25, 0.25]

################# Fitting Model, Saving Simulation Results of Interest #################

M = mc.MCMC(locals())

# 20000 draws, only consider the final 2000 draws for determining cluster assignment though
M.sample(20000, 0)


# For now we are only saving the last 2000 draws of the category parameter
# The trace for other parameters may be obtained and saved as well 
categorydraw = M.trace('category').gettrace()

np.savetxt('cat_4var_fit3.csv', categorydraw, delimiter=',')


