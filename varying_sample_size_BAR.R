##### File #####################################################################
# Author: Luke Ondijo Ouma
# Last edited: 0/04/2020
#
# Details: File to evaluate operational characteristics for BAR approaches of Tx 
#          allocation for patients with multiple biomarkers.
#         - Here BAR is used to guide the allocation of patients with Multiple biomarkers OR 
#           eligible for multiple subgroups

################# load relevant packages ##########
if(!require("rjags")) install.packages("rjags"); library(rjags)
if(!require("R2jags")) install.packages("R2jags"); library(R2jags)

if(!require("BBmisc")) install.packages("BBmisc"); library(BBmisc)
#if(!require("devtools")) install.packages("devtools"); library(devtools)
if(!require("coda")) install.packages("coda"); library(coda)
if(!require("R2OpenBUGS")) install.packages("R2OpenBUGS"); library(R2OpenBUGS)
if(!require("boot")) install.packages("boot"); library(boot)

# devtools::install_version("elrm", "1.2.4") 
library(elrm)
# if(!require(foreach)) install.packages("foreach"); library(foreach)
# if(!require(doParallel)) install.packages("doParallel"); library(doParallel)
if(!require(openxlsx)) install.packages("openxlsx"); library(openxlsx)
###################################################
setwd("C:/Users/Luke Ondijo/OneDrive - Newcastle University/PhD research work/Multiple Biomarkers project/4_BAR/JAGS_1-6")

#------------------------------------------------------------------------------------------------------#
#
#     Variable initialisation   
# npatients                  <- 400                     # No of patients
# prevalence                 <- c(0.3, 0.25, 0.3, 0.25) # Biomarker prevalence
#parameters:
# K - number of treatments
# biomarkerprevalences - K dimensional vector of prevalence of each biomarker. Biomarker positivity is assumed to 
# n.eachstage - vector of length J (number of stages), with the timing of interim analyses, in terms of recruitment
# alpha       - intercept in response simulation model
# beta        - K-dimensional vector of treatment effects (on lor scale)
# gamma       - K-dimension vector of marginal biomarker effects
# delta       - KxK matrix of biomarker-treatment interactions; i,jth entry represents effect when patient positive for biomarker j is given experimental treatment i
# a and b     - tuning parameters for adaptive randomization procedure
# maxN        - total number of patients
# delay       - delay in months between recruitment and assessment
# recruitmentpermonth - number of patients recruited per month on average

rowMax <- function(mat)
{
  return(apply(mat,1,max))
}


getposteriormeans_jags <- function(biomarkerprofiles,treatmentvector,responses,allprofiles,K,N, datafilename)
{
  # r         <- responses
  # biomarker <- biomarkerprofiles
  # # N         <- length(tempallocation)
  # N         <- N
  # K         <- K
  # allprofiles = allprofiles
  
  #define covariates representing treatment:
  treatmentmatrix <- matrix(0,N, K)
  for(k in 1:K)
  {
    treatmentmatrix[,k] <- ifelse(treatmentvector==k+1,1,0)
  }
  
  datastring=paste("list( biomarker=structure(.Data=c(",collapse(as.double(t(biomarkerprofiles))),"),.Dim=c(",N,",",K,")),treatment=structure(.Data=c(",collapse(as.double(t(treatmentmatrix))),"),.Dim=c(",N,",",K,")),allprofiles=structure(.Data=c(",collapse(as.double(t(allprofiles))),"),.Dim=c(",length(allprofiles[,1]),",",K,")),r=c(",collapse(as.double(responses)),"),K=",K,",N=",N,")",sep="")
  
  cat(datastring,file=datafilename)
  
  # bugs.data(datafilename, data.file = paste(datafilename, ".txt", sep = ""))
  # Inits
  inits <- function(){
    list(
      alpha= 0, beta= c(0,0,0,0),gamma=c(0,0,0,0),delta=matrix(0, 4,4)) 
  }
  
  parameters = c("postprob_betterthancontrol")
  
  MCMCSim1 <- jags(data =dget(datafilename), inits = inits, parameters.to.save= parameters, model.file = "model_K=4_betterthancontrol.txt", 
                   n.chains = 1, n.iter = 5000, n.burnin = 2500, progress.bar = "none", working.directory = "C:/Users/Luke Ondijo/OneDrive - Newcastle University/PhD research work/Multiple Biomarkers project/4_BAR/JAGS_1-6")
  # dat = dget("data_2.97873857896775.txt")
  
  #return as matrix with each row corresponding to a row of allprofiles
  
  # return(list(allocation_betterthancontrol = MCMCSim1$mean$postprob_betterthancontrol))
  return(list(allocation_betterthancontrol = t(matrix(MCMCSim1$BUGSoutput$mean$postprob_betterthancontrol, K, length(allprofiles[,1])))))
  
}
#                          FUNCTION         
#-------------------------------------------------------------------------------------------------------#
simulation_replicate_jags <- function(alpha,beta,gamma,delta,biomarkerprevalences,K,a,b,maxN,n.eachstage,delay,recruitmentpermonth)
{
  
  # nbiomarkers                <- length(biomarkerprevalences)
  # ntreatments                <- nbiomarkers + 1L    # Control and 4 experimental treatments
  # K                          <- nbiomarkers         # No experimental treatments
  
  #set number of stages:    
  J <- length(n.eachstage)
  
  # Find number of patients to be recruited between each analysis
  recruiteachstage <-  c(n.eachstage[1],n.eachstage[-1]-n.eachstage[-length(n.eachstage)])
  
  datafilename     <- paste("data_",runif(1,0,10), ".txt",sep="") #####What does this actually do?
  
  # Determine the times of interim analysis
  incrementrecruitmenttimes <- rexp(maxN, recruitmentpermonth)
  recruitmenttimes          <- cumsum(incrementrecruitmenttimes)
  # Time of interim analysis
  interimanalysistimes      <- recruitmenttimes[n.eachstage]
  # add small amount on to initial J-1 analysis times
  interimanalysistimes[1:(J-1)] <- interimanalysistimes[1:(J-1)] + 0.0001
  interimanalysistimes[J]       <- interimanalysistimes[J] + delay + 0.0001
  
  #----------- Obtain patient biomarker profiles
  seed= as.integer(runif(1,1,1e4)*1.856 +169254)
  set.seed(seed)
  
  allprofiles <- NULL
  
  for(i1 in 0:1){    for(i2 in 0:1){      for(i3 in 0:1){        {for (i4 in 0:1) {
    {allprofiles <- rbind(allprofiles,c(i4,i3,i2,i1))}}}}}}
  
  
  # Determine all patients biomarker profiles:
  biomarkerprofiles <- matrix(0, maxN, K)
  
  for(i in 1:maxN)
  {
    biomarkerprofiles[i,] <- rbinom(K, 1, biomarkerprevalences)
  }
  
  # Convert to binary score:
  binaryscore   <- rep(1, maxN)
  for(i in 1:K)
  {
    binaryscore <- binaryscore + biomarkerprofiles[,i]*2^{i-1}
  }
  
  # Get patient numbers recruited at each stage:
  patientranges        <- matrix(0, J, 2)
  patientranges[1,]    <- c(1, n.eachstage[1])
  for(i in 2:J)
  {
    patientranges[i,]  <- c(n.eachstage[i-1]+1,n.eachstage[i])
  }
  
  # Set initial allocation probabilities
  allocationprobabilities          <- matrix(0, 2^K, K+1)
  allocationprobabilities[1,]=rep(1/5,5)
  allocationprobabilities[2,]=c(0.5,0.5,0,0,0)
  allocationprobabilities[3,]=c(0.5,0,0.5,0,0)
  allocationprobabilities[4,]=c(1/3,1/3,1/3,0,0)
  allocationprobabilities[5,]=c(1/2,0,0,1/2,0)
  allocationprobabilities[6,]=c(1/3,1/3,0,1/3,0)
  allocationprobabilities[7,]=c(1/3,0,1/3,1/3,0)
  allocationprobabilities[8,]=c(1/4,1/4,1/4,1/4,0)
  allocationprobabilities[9,]=c(0.5,0,0,0,0.5)
  allocationprobabilities[10,]=c(1/3,1/3,1/3,0,0)
  allocationprobabilities[11,]=c(1/3,0,1/3,0,1/3)
  allocationprobabilities[12,]=c(1/4,1/4,1/4,0,1/4)
  allocationprobabilities[13,]=c(1/3,0,0,1/3,1/3)
  allocationprobabilities[14,]=c(1/4,1/4,0,1/4,1/4)
  allocationprobabilities[15,]=c(1/4,0,1/4,1/4,1/4)
  allocationprobabilities[16,]=c(1/5,1/5,1/5,1/5,1/5)
  #allocationprobabilities         <- matrix(1/(K+1), 2^K, K+1)
  allocationprobabilities_init    <- allocationprobabilities
  posteriormeans_betterthancontrol <- allocationprobabilities[,-1]
  
  allocation              <- rep(0, maxN)
  response                <- rep(0, maxN)
  allocationpersubgroup   <- matrix(0, 2^K, K+1)
  diff                    <- rep(0, 2^K)
  patientsonbesttreatment <- numeric(length = J)
  
  for (n.stage in 1:J) {
    
    # Get allocation and response vectors for this stages patients
    tempallocation        <- rep(0, recruiteachstage[n.stage])
    tempresponse          <- rep(0, recruiteachstage[n.stage])
    
    # Get binary score for each patient
    tempbinaryscore       <- binaryscore[patientranges[n.stage,1]:patientranges[n.stage,2]]
    tempbiomarkerprofiles <- biomarkerprofiles[patientranges[n.stage,1]:patientranges[n.stage,2],]
    tempmbiomarkerpts     <- which(rowSums(tempbiomarkerprofiles)>1)
    
    if(n.stage == 1)
    {
      for (i in 1:length(tempbinaryscore)) {
        
        tempallocation[i]  <- which(as.double(rmultinom(1,1,allocationprobabilities[tempbinaryscore[i],]))==1)
        allocationpersubgroup[tempbinaryscore[i], tempallocation[i]] <- allocationpersubgroup[tempbinaryscore[i], tempallocation[i]] + 1
        diff               <- rowMax(allocationpersubgroup[,(2:4)]) - allocationpersubgroup[,1]
      }
    }
    
    if(n.stage >1)
    {
      for (i in 1:length(tempbinaryscore)) {
        
        allocationprobabilities.temp     <- cbind(rep(0,length(allprofiles[,1])), posteriormeans_betterthancontrol^ (a*((patientranges[n.stage,1]+i)/maxN)^b))
        allocationprobabilities.temp[,1] <- apply(allocationprobabilities.temp, 1, max)* (exp(diff)^ ((1/(K+1))* (patientranges[n.stage,1] +i)/maxN))
        allocationprobabilities          <- allocationprobabilities.temp/rowSums(allocationprobabilities.temp)
        
        if(i %in% tempmbiomarkerpts){   # Patients with multiple biomarkers, the probabilities are updated here
          tempallocation[i]  <- which(as.double(rmultinom(1,1,allocationprobabilities[tempbinaryscore[i],]))==1)
          #allocationpersubgroup[tempbinaryscore[i], tempallocation[i]] <- allocationpersubgroup[tempbinaryscore[i], tempallocation[i]] + 1
          #diff               <- rowMax(allocationpersubgroup[,(2:4)]) - allocationpersubgroup[,1] 
        }else {           # # Patients without multiple biomarkers, the probabilities are just same a sthe initial
          tempallocation[i]  <- which(as.double(rmultinom(1,1,allocationprobabilities_init[tempbinaryscore[i],]))==1)
          #allocationpersubgroup[tempbinaryscore[i], tempallocation[i]] <- allocationpersubgroup[tempbinaryscore[i], tempallocation[i]] + 1
          #diff               <- rowMax(allocationpersubgroup[,(2:4)]) - allocationpersubgroup[,1]
        }
        
        # tempallocation[i]  <- which(as.double(rmultinom(1,1,allocationprobabilities[tempbinaryscore[i],]))==1)
        allocationpersubgroup[tempbinaryscore[i], tempallocation[i]] <- allocationpersubgroup[tempbinaryscore[i], tempallocation[i]] + 1
        diff               <- rowMax(allocationpersubgroup[,(2:4)]) - allocationpersubgroup[,1]
        
      }
    }
    
    # Obtain logodds for each patient
    logodds.response <- sapply(1:length(tempallocation), function(x){
      return(alpha + ifelse(tempallocation[x]>1, beta[tempallocation[x]-1] + as.double(delta[tempallocation[x]-1,] %*% tempbiomarkerprofiles[x,]),0) + tempbiomarkerprofiles[x,] %*% gamma)
      
    })
    
    #patientsonbesttreatment[n.stage] = if(max(logodds.response) != min(logodds.response)){length(which(inv.logit(logodds.response)== max(inv.logit(logodds.response))))}
    # if(max(logodds.response) != min(logodds.response)){patientsonbesttreatment[n.stage] = length(which(inv.logit(logodds.response)== max(inv.logit(logodds.response))))}else {
    #   patientsonbesttreatment[n.stage] =0
    # }
    patientsonbesttreatment[n.stage] =ifelse(max(logodds.response) != min(logodds.response),length(which(inv.logit(logodds.response)== max(inv.logit(logodds.response)))),0)
    # Simulate response/Outcome
    # First obtain probability of response
    tempresponse <- rbinom(length(tempallocation),1, inv.logit(logodds.response))
    
    # Fill in allocation and response
    response[patientranges[n.stage,1]:patientranges[n.stage,2]]   <- tempresponse
    allocation[patientranges[n.stage,1]:patientranges[n.stage,2]] <- tempallocation
    
    # Fit Bayesian (BUGS) model to get new allocation probabilities
    # First obtain data for patients who are assessed by time of interim analysis
    if(n.stage < J)
    {
      tempbiomarkerprofiles <- biomarkerprofiles[(recruitmenttimes + delay) < interimanalysistimes[n.stage],]
      tempallocation        <- allocation[(recruitmenttimes + delay) < interimanalysistimes[n.stage]]
      tempresponse          <- response[(recruitmenttimes + delay) < interimanalysistimes[n.stage]]
      
      posteriormeans                   <-  getposteriormeans_jags(tempbiomarkerprofiles, tempallocation, tempresponse, allprofiles, K, length(tempallocation), datafilename)
      
      posteriormeans_betterthancontrol <- posteriormeans$allocation_betterthancontrol + 0.00001 #Add tiny amount to avoid problems with all 0s
      
    }
    
  }
  
  #-------- Final analysis
  
  glm1 <- glm(response ~ as.factor(allocation)*biomarkerprofiles, family = "binomial")
  
  # Get variance for testing first experimental treatment in patients with first biomarker positive:
  # Get matrix of test statistics:
  temp <- matrix(0, 25, 20)
  #first get tests of treatment 1 in different groups:
  temp[,1] = c(0, 1, rep(0,23))
  temp[,2] = c(0, 1, rep(0,7), 1, rep(0, 15))
  temp[,3] = c(0, 1, rep(0,11), 1, rep(0, 11))
  temp[,4] = c(0, 1, rep(0,15), 1, rep(0, 7))
  temp[,5] = c(0, 1, rep(0,19), 1, rep(0, 3))
  
  # tests of treatment 2 in different groups:
  temp[,6] = c(rep(0,2), 1, rep(0,22))
  temp[,7] = c(rep(0,2), 1, rep(0,7), 1, rep(0, 14))
  temp[,8] = c(rep(0,2), 1, rep(0,11), 1, rep(0, 10))
  temp[,9] = c(rep(0,2), 1, rep(0,15), 1, rep(0, 6))
  temp[,10] = c(rep(0,2), 1, rep(0,19), 1, rep(0, 2))
  
  # tests of treatment 3 in different groups:
  temp[,11] = c(rep(0,3), 1, rep(0,21))
  temp[,12] = c(rep(0,3), 1, rep(0,7), 1, rep(0, 13))
  temp[,13] = c(rep(0,3), 1, rep(0,11), 1, rep(0, 9))
  temp[,14] = c(rep(0,3), 1, rep(0,15), 1, rep(0, 5))
  temp[,15] = c(rep(0,3), 1, rep(0,19), 1, 0)
  
  # tests of treatment 4 in different groups:
  temp[,16] = c(rep(0,4), 1, rep(0,20))
  temp[,17] = c(rep(0,4), 1, rep(0,7), 1, rep(0, 12))
  temp[,18] = c(rep(0,4), 1, rep(0,11), 1, rep(0, 8))
  temp[,19] = c(rep(0,4), 1, rep(0,15), 1, rep(0, 4))
  temp[,20] = c(rep(0,4), 1, rep(0,19), 1)
  
  # Check if any coefficients are NA
  whichna   <- is.na(as.double(glm1$coef))
  varmatrix <- summary(glm1)$cov.unscaled
  
  if(any(whichna))
  {
    # set corresponding coefficient to 0
    glm1$coef[whichna]             <- 0
    
    # Add a row and column to glm1$cov.unscaled with infinite entries
    tempmatrix <- matrix(0, 25, 25)
    tempmatrix[!whichna, !whichna] <- summary(glm1)$cov.unscaled
    tempmatrix[whichna,]           <- 9e50
    tempmatrix[, whichna]          <- 9e50
    varmatrix                      <- tempmatrix
  }
  
  mean            <- as.double(t(temp) %*% glm1$coef)
  var             <- t(temp) %*% varmatrix %*% temp
  teststatistics  <- mean/sqrt(diag(var))
  onesidedpvalues <- pnorm(-abs(teststatistics))
  
  
  # Number/Proportion of patients being treated
  patientsontreatment <- length(which(allocation !=1))/ maxN *100
  # Number allocated to each arm
  patients_per_arm   <- table(allocation)
  #proportionpatients_per_arm   <- prop.table(patients_per_arm)
  
  # Number responded on each arm
  responders_per_arm <- table(response, allocation)
  #proportionresponders_per_arm <- prop.table(responders_per_arm, margin = 2)
  
  # Number of patients on the best treatment available to them (patientsonbesttreatment)
  #- This maybe interesting because we have to consider stagewise
  #- Also a patient maybe on the best treatment twice, at stage 1 and 2 or 1 and 4
  # - May want to think of how to acocunt for that.
  # Alternative - how many pts are on the best treatment on average at each stage?
  #- we consider the best treatment as the one that maximizes the probability of response for a patient
  patientsonBesttreatment = sum(patientsonbesttreatment)/maxN*100
  
  # Delete data file
  system(paste("rm ", datafilename, sep = ""))
  
  #Return results
  return(list(teststatistics=teststatistics, onesidedpvalues= onesidedpvalues, allocation=allocation, response=response, biomarkerprofiles=biomarkerprofiles, 
              allocationpersubgroup= allocationpersubgroup, patientsontreatment= patientsontreatment, 
              patients_per_arm= patients_per_arm, responders_per_arm= responders_per_arm, patientsonBesttreatment= patientsonBesttreatment))
  
  # end code 
}

#Load simulation scenarios
source("C:/Users/Luke Ondijo/OneDrive - Newcastle University/PhD research work/Multiple Biomarkers project/simulation_scenarios_BAR.R")

biomarkerprevalences  <- c(0.3, 0.25, 0.3, 0.25)
K                     <- 4
#n.eachstage           <- nn.eachstage = c(100,175,250,325,400) 
recruitmentpermonth   <- 8
delay                 <- 6
#maxN                  <- 400 maxN= samplesizes[1]
niterations           <- 5        # 3000
a=13.5; b=2.75
#maxN= samplesizes[1]; n.eachstage= numeachstage[[1]]
#------------------- Run the code

samplesizes =seq(100,1000, by=50)
numeachstage = list()
for(i in 1:length(samplesizes)){
  stage1 = round(0.25*samplesizes[i],0)
  next_stage = round(0.25*(0.75*samplesizes[i]),0)
  numeachstage[[i]] <- c(stage1, stage1+next_stage, stage1+ 2*next_stage, stage1+3*next_stage,samplesizes[i])
}


library(parallel)
c1     <- makeCluster(detectCores())
clusterExport(c1, varlist = c("samplesizes","sim_scenario", "getposteriormeans_jags", "simulation_replicate_jags", "rowMax",
                              "biomarkerprevalences","K","a", "b","numeachstage","delay","recruitmentpermonth","niterations"))
# clusterExport(c1, varlist = c("sim_scenario", 'getposteriormeans_jags', 'simulation_replicate_jags', 'rowMax',
#                               'biomarkerprevalences','K','a','b','maxN','n.eachstage','delay','recruitmentpermonth', "niterations"))
clusterEvalQ(c1, library(parallel))
#clusterEvalQ(c1, library(doParallel))
clusterEvalQ(c1, library(boot))
clusterEvalQ(c1, library(data.table))
clusterEvalQ(c1, library(Rfast))
clusterEvalQ(c1, library(dplyr))
clusterEvalQ(c1, library(rjags))
clusterEvalQ(c1, library(R2jags))
clusterEvalQ(c1, library(elrm))
clusterEvalQ(c1, library(BBmisc))

# start = Sys.time()
# BAR_results_JAGS_15      <-  parLapply(c1, 1:5, function(i){    # Evaluate 6 scenarios 
#   replicate(niterations, 
#             simulation_replicate_jags(alpha = sim_scenario[[i]]$alpha,  beta = sim_scenario[[i]]$beta, gamma = sim_scenario[[i]]$gamma,  delta= sim_scenario[[i]]$delta,
#                                       biomarkerprevalences,K,a,b,maxN,n.eachstage,delay,recruitmentpermonth), simplify = FALSE)
# })
# stopCluster(c1)
# 
# stop = Sys.time()
# print(paste("Run time is:", stop-start,sep=""))

start = Sys.time()
# BAR_results_JAGS_varyssize_test2      <-  try(parLapply(c1, 4:7, function(k) {    #length(samplesizes)
# #  lapply(c(1:2,8), function(i) {
#   replicate(niterations,
#             simulation_replicate_jags(alpha = sim_scenario[[1]]$alpha,  beta = sim_scenario[[1]]$beta, gamma = sim_scenario[[1]]$gamma,  delta= sim_scenario[[1]]$delta,
#                                       biomarkerprevalences=biomarkerprevalences, K=K, a=a, b=b, maxN =samplesizes[k],
#                                       n.eachstage= numeachstage[[k]], delay=delay, recruitmentpermonth=recruitmentpermonth), simplify = FALSE)
# }))
# #})
start = Sys.time()
BAR_results_JAGS_varyssize_test = replicate(niterations,
       simulation_replicate_jags(alpha = sim_scenario[[1]]$alpha,  beta = sim_scenario[[1]]$beta, gamma = sim_scenario[[1]]$gamma,  delta= sim_scenario[[1]]$delta,
       biomarkerprevalences=biomarkerprevalences, K=K, a=a, b=b, maxN =samplesizes[4],
       n.eachstage= numeachstage[[4]], delay=delay, recruitmentpermonth=recruitmentpermonth), simplify = FALSE)

#stopCluster(c1)
stop = Sys.time()
print(paste("Run time is:", stop-start,sep=""))

# Save results as RDS file to a folder
saveRDS(BAR_results_JAGS_varyssize, "C:/Users/Luke Ondijo/OneDrive - Newcastle University/PhD research work/Multiple Biomarkers project/results/BAR_results_JAGS_varyssize")
