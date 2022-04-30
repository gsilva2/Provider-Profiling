library(cluster)
library(ggplot2)
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

############## Clustering Analyses ##############

########### Obtaining clusters for each value of M ###########

# Reading in cluster posterior draws obtained from running script2_simdata_fitn_2cont2cat.py
# These full scripts have not been included in GitHub because they are very large files 
# The cluster assignments that result from running this section for M=2, 3, 4, 5 has been included ("cluster_assignment_M.Rdata")
#dat_2<-as.matrix(read.csv("cat_4var_fit2.csv", header=FALSE))
#dat_3<-as.matrix(read.csv("cat_4var_fit3.csv", header=FALSE))
#dat_4<-as.matrix(read.csv("cat_4var_fit4.csv", header=FALSE))
#dat_5<-as.matrix(read.csv("cat_4var_fit5.csv", header=FALSE))

#iter<-18000:20000

#clus_2<-c()
#clus_3<-c()
#clus_4<-c()
#clus_5<-c()
#for(i in 1:dim(dat_3)[2]){
#  clus_2<-c(clus_2, as.numeric(names(which.max(table(dat_2[iter,i]))))+1)
#  clus_3<-c(clus_3, as.numeric(names(which.max(table(dat_3[,i]))))+1)
#  clus_4<-c(clus_4, as.numeric(names(which.max(table(dat_4[iter,i]))))+1)
#  clus_5<-c(clus_5, as.numeric(names(which.max(table(dat_5[iter,i]))))+1)
#}

#clus_data<-as.data.frame(cbind(prov=1:n_prov, clus_2, clus_3, clus_4, clus_5))
#save(clus_data, file="cluster_assignment_M.Rdata")
load(file="cluster_assignment_M.Rdata")

########### Identifying M* ###########

# Generating distance matrix
sum_num<-function(dataset, variable){
  n<-which(names(dataset)==variable)
  one<-tapply(dataset[,n], dataset$Provider, mean)
  return(one)
}

final<-as.data.frame(cbind(sum_num(data, "age"), sum_num(data, "hosp_los")))

final<-cbind(final, sum_num(data, "female"), sum_num(data, "renal"))

# Rescaling so that all entries in the distance matrix are weighted equally
for(j in 1:dim(final)[2]){
  final[,j]<-scales::rescale(final[,j], to=c(0, 1))
}

# Computing distance matrix
X<-as.matrix(final)
X_dist<-as.matrix(dist(X, method = "euclidean"))

hi<-silhouette(clus_data$clus_5, dmatrix=X_dist)
avg_5<-mean(hi[ ,3])
hi<-silhouette(clus_data$clus_4, dmatrix=X_dist)
avg_4<-mean(hi[ ,3])
hi<-silhouette(clus_data$clus_3, dmatrix=X_dist)
avg_3<-mean(hi[,3])
hi<-silhouette(clus_data$clus_2, dmatrix=X_dist)
avg_2<-mean(hi[,3])

clus_plot_dat<-as.data.frame(cbind(x=seq(2, 5), y=c(avg_2, avg_3, avg_4, avg_5)))
g<-ggplot(data =clus_plot_dat, mapping=aes(x=x, y=y))+geom_line(size=1)+geom_point(size=4)
g<-g+labs(x="Number of Clusters", y="Average Silhouette Width")
g

# From this plot we see that the value of M that maximizes the average silhouette width is 3
# M*=3 for this example
