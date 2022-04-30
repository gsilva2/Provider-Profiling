library(MASS)
library(boot)
library(lme4)

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

############## Outcome Data Analyses ##############

########### Merging clustering information with dataset ###########

# No clustering grouping (NoCLUS)
# All providers belong to the same cluster
NoClus<-rep(1, times=dim(data)[1])

# Actual clustering grouping (CLUS)
# This is obtained obtained after clustering the providers in script_simdata_fitn_2cont2cat.py and using the results on M* and the clusters each provider had been assigned to
# From script3_clusteranalysis.R, M*=3 and so we use this column for clustering
#load(file="cluster_assignment_M.Rdata")
clus<-clus_data$clus_3
hi<-t(rbind(1:n_prov, clus))
colnames(hi)<-c("Provider", "Clus")
data<-merge(data, hi, by="Provider")
Clus<-data$Clus

# Random clustering grouping (RanCLUS)
# This is random but was done in this manner to make sure that each of the random groups has the same number of providers as clustering and that each cluster has 6 outliers to identify
ranclus<-rep(c(1, 2, 3), length=120)
hi2<-t(rbind(1:n_prov, ranclus))
colnames(hi2)<-c("Provider", "RanClus")
data<-merge(data, hi2, by="Provider")
RanClus<-data$RanClus

########### Estimation ###########

###### Function to estimate SRRs ###### 
# clus_group is the clustering vector that describes if all the providers are in the same group (NoCLUS) or if they are in different groups (CLUS, RanCLUS)
srr<-function(dat, clus_group){
  
  res_final<-list()
  for(m in 1:length(dat)){
    
    work_dat<-dat[[m]]
    work_dat<-as.data.frame(cbind(work_dat, clus_group))
    
    res_mat<-list()
    for(k in 1:length(table(clus_group))){
      #Fitting the model within the cluster
      data_fit<-work_dat[work_dat$clus_group==k, ]
      data_fit$Provider<-factor(data_fit$Provider)
      mod<-glmer(sim_Y~age+hosp_los+female+renal+(1|Provider), family="binomial", data=data_fit)
      
      #### Splitting data for when we estimate SRRs
      sub<-split(data_fit, data_fit$Provider)  
      
      res<-matrix(NA, nrow=length(sub), ncol=3)
      for(i in 1:length(sub)){
        data_est<-sub[[i]]
        #Putting data in proper matrix form for predictions
        data_est<-data_est[ ,2:5]
        data_est<-as.matrix(cbind(rep(1, times=dim(data_est)[1]), data_est))
        
        intercept<-coef(mod)$Provider[i, ][[1]]
        coefficients<-as.matrix(coef(mod)$Provider[i, ])
        pred<-inv.logit(data_est%*%t(coefficients))
        obs<-sum(pred)
        coefficients<-as.matrix(summary(mod)$coefficients[,1])
        pred<-inv.logit(data_est%*%coefficients)
        exp<-sum(pred)
        estim<-obs/exp
        res[i, ]<-c(as.numeric(as.character(sub[[i]]$Provider[1])), estim, intercept)
      }
      res_mat[[k]]<-res
    }
    
    res1<-c()
    res2<-c()
    res3<-c()
    for(r in 1:length(table(clus_group))){
      res1<-c(res1, res_mat[[r]][,1])
      res2<-c(res2, res_mat[[r]][,2])
      res3<-c(res3, res_mat[[r]][,3])
    }
    finale<-data.frame(res1, res2, res3, noclus=rep(1, n_prov), clus, ranclus)
    finale<-finale[order(res1), ]
   
    res_final[[m]]<-finale
    
    print(m)
  }
  
  return(res_final)
  
}


###### Estimating under different groupings, configurations ###### 
# Here only logprop, logimp2 specification is considered but other simulation settings could also be explored
# load(file="logprop_dat.Rdata")
# load(file="logimp2_dat.Rdata")
logprop_1clus<-srr(log_prop, NoClus)
logprop_clus<-srr(log_prop, Clus)
logprop_clusrand<-srr(log_prop, RanClus)

logimp1_1clus<-srr(log_imp1, NoClus)
logimp1_clus<-srr(log_imp1, Clus)
logimp1_clusrand<-srr(log_imp1, RanClus)

########### Analysis ###########
#load(file="outliers_sim_15.RData")
out_ind<-which(outliers==TRUE)
notout_ind<-setdiff(as.numeric(1:n_prov), out_ind)

###### Function to extract FPR, FNR, Accuracy ###### 
dat_ext<-function(list_sav, clus_type){
  m<-length(list_sav)
  res_mat<-matrix(NA, nrow=m, ncol=3)
  
  fnr<-c()
  fpr<-c()
  acc<-c()
  
  for(i in 1:m){
    #These need to be modified depending on number of outliers
    da<-list_sav[[i]]
    
      l<-which(names(da)==clus_type)
      new_dat<-data.frame(da[,1], da[,2], da[,l])
      names(new_dat)<-c("prov", "srr", "clus")
      new_dat$clus<-as.factor(new_dat$clus)
      
      out_method<-c()
      for(j in 1:length(table(new_dat$clus))){
        num_in_clus<-table(new_dat$clus)[j]
        #The percent of outliers in each cluster was assumed to be 15% so let's calculate the worst 15% within each cluster
        num_out<-0.15*num_in_clus
        min_dat<-new_dat[new_dat$clus==j, ]
        min_rank<-rank(min_dat[,2])
        out_method<-c(out_method, min_dat[min_rank>num_in_clus-num_out, 1])
        }
      
      
      notout_method<-setdiff(1:n_prov, out_method)
      
      out_method<-1:n_prov%in%out_method
      notout_method<-1:n_prov%in%notout_method
      
      fnr<-c(fnr, (1-(sum(out_method&outliers)/sum(outliers)))*100)
      fpr<-c(fpr, (1-(sum(notout_method&!outliers)/sum(!outliers)))*100)
      acc<-c(acc, ((sum(out_method&outliers)+sum(notout_method&!outliers))/n_prov)*100)
  }
  res_mat[, 1]<-fnr
  res_mat[, 2]<-fpr
  res_mat[, 3]<-acc
  return(res_mat)
}

###### Function to summarize FPR, FNR, Accuracy ###### 
sum_dat<-function(dat_extract){
  med<-round(apply(dat_extract, 2, median), 2)
  ci1<-round(apply(dat_extract, 2, quantile, 0.025), 2)
  ci2<-round(apply(dat_extract, 2, quantile, 0.975), 2)
  res<-paste(med, " (", ci1, ", ", ci2, ")", sep="")
  return(res)
}


#modify name of dataset accordingly to correspond to res1, res2, res3, res4, res5, res6
one<-sum_dat(dat_ext(logprop_1clus, "noclus"))
two<-sum_dat(dat_ext(logprop_clus, "clus"))
three<-sum_dat(dat_ext(logprop_clusrand, "ranclus"))

res1<-c(one[1], two[1], three[1], one[2], two[2], three[2], one[3], two[3], three[3])

one<-sum_dat(dat_ext(logimp1_1clus, "noclus"))
two<-sum_dat(dat_ext(logimp1_clus, "clus"))
three<-sum_dat(dat_ext(logimp1_clusrand, "ranclus"))

res2<-c(one[1], two[1], three[1], one[2], two[2], three[2], one[3], two[3], three[3])
res<-as.data.frame(cbind(res1, res2))
rownames(res)<-c("FNR NoCLUS", "FNR CLUS", "FNR RanCLUS", "FPR NoCLUS", "FPR CLUS", "FPR RanCLUS", "Accuracy NoCLUS", "Accuracy CLUS", "Accuracy RanCLUS")
colnames(res)<-c("Log-Prop", "Log-Imp1")
