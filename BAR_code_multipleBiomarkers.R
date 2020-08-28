##### File #####################################################################
# Author: Luke Ouma
# Last edited: 0/04/2020
#
# Details: File to evaluate operational characteristics for BAR approache of Tx allocation for patients with multiple biomarkers.
#         - Here BAR is used to guide the allocation of patients with Multiple biomarkers only

################# load relevant packages ##########
if(!require("rjags")) install.packages("rjags"); library(rjags)
if(!require("BBmisc")) install.packages("BBmisc"); library(BBmisc)
if(!require("devtools")) install.packages("devtools"); library(devtools)
if(!require("coda")) install.packages("coda"); library(coda)
if(!require("R2OpenBUGS")) install.packages("R2OpenBUGS"); library(R2OpenBUGS)
if(!require("boot")) install.packages("boot"); library(boot)

devtools::install_version("elrm", "1.2.4") 
library(elrm)
if(!require(openxlsx)) install.packages("openxlsx"); library(openxlsx)

###################################################
setwd("Multiple Biomarkers project/4_BAR/BAR_712")

#------------------------------------------------------------------------------------------------------#
#
#     Variable initialisation   
# npatients   - No of patients
# prevalence  - Biomarker prevalence

# Parameters:
# K           - number of treatments
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

  #define covariates representing treatment:
  treatmentmatrix <- matrix(0,N, K)
  for(k in 1:K)
  {
    treatmentmatrix[,k] <- ifelse(treatmentvector==k+1,1,0)
  }
  
  datastring=paste("list( biomarker=structure(.Data=c(",collapse(as.double(t(biomarkerprofiles))),"),.Dim=c(",N,",",K,")),treatment=structure(.Data=c(",collapse(as.double(t(treatmentmatrix))),"),.Dim=c(",N,",",K,")),allprofiles=structure(.Data=c(",collapse(as.double(t(allprofiles))),"),.Dim=c(",length(allprofiles[,1]),",",K,")),r=c(",collapse(as.double(responses)),"),K=",K,",N=",N,")",sep="")
  
  cat(datastring,file=datafilename)
  
  # Inits
  inits <- function(){
    list(
      alpha= 0, beta= c(0,0,0,0),gamma=c(0,0,0,0),delta=matrix(0, 4,4)) 
  }
  
  parameters = c("postprob_betterthancontrol")
  
  MCMCSim1 <- jags(data =dget(datafilename), inits = inits, parameters.to.save= parameters, model.file = "model_K=4_betterthancontrol.txt", 
                   n.chains = 1, n.iter = 5000, n.burnin = 2500, progress.bar = "none", working.directory = "Multiple Biomarkers project/4_BAR/BAR_712")
  
    return(list(allocation_betterthancontrol = t(matrix(MCMCSim1$BUGSoutput$mean$postprob_betterthancontrol, K, length(allprofiles[,1])))))
  
}
#                          FUNCTION         
#-------------------------------------------------------------------------------------------------------#
simulation_replicate_jags <- function(alpha,beta,gamma,delta,biomarkerprevalences,K,a,b,maxN,n.eachstage,delay,recruitmentpermonth)
{
  
  #set number of stages:    
  J <- length(n.eachstage)
  
  # Find number of patients to be recruited between each analysis
  recruiteachstage <-  c(n.eachstage[1],n.eachstage[-1]-n.eachstage[-length(n.eachstage)])
  
  datafilename     <- paste("data_",runif(1,0,10), ".txt",sep="") 
  
  # Determine the times of interim analysis
  incrementrecruitmenttimes <- rexp(maxN, recruitmentpermonth)
  recruitmenttimes          <- cumsum(incrementrecruitmenttimes)
  
  # Time of interim analysis
  interimanalysistimes      <- recruitmenttimes[n.eachstage]
  
  # add small amount on to initial J-1 analysis times
  interimanalysistimes[1:(J-1)] <- interimanalysistimes[1:(J-1)] + 0.0001
  interimanalysistimes[J]       <- interimanalysistimes[J] + delay + 0.0001
  
  # Obtain patient biomarker profiles
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

        }else {           # # Patients without multiple biomarkers, the probabilities are just same a sthe initial
          tempallocation[i]  <- which(as.double(rmultinom(1,1,allocationprobabilities_init[tempbinaryscore[i],]))==1)

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
    

    patientsonbesttreatment[n.stage] =ifelse(max(logodds.response) != min(logodds.response),length(which(inv.logit(logodds.response)== max(inv.logit(logodds.response)))),0)
    
    # Simulate response/Outcome
    ## First obtain probability of response
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

# Load simulation scenarios
source("Multiple Biomarkers project/simulation_scenarios_BAR.R")

biomarkerprevalences  <- c(0.3, 0.25, 0.3, 0.25)
K                     <- 4
n.eachstage           <- c(100,175,250,325,400)
recruitmentpermonth   <- 8
delay                 <- 6
maxN                  <- 400
niterations           <- 10000        # 2500
a=13.5; b=2.75

# Execute sim in parallel
library(parallel)
c1     <- makeCluster(detectCores())

clusterExport(c1, varlist = c("sim_scenario", 'getposteriormeans_jags', 'simulation_replicate_jags', 'rowMax',
                              'biomarkerprevalences','K','a','b','maxN','n.eachstage','delay','recruitmentpermonth', "niterations"))
clusterExport(c1, varlist = c('multiple_biomarker_allocation', 'constrained_randomization', 
                              'equal_randomization_eligible_Tx', 'hierarchy', 'prioritize_randomization_eligible_Tx'))
clusterEvalQ(c1, library(parallel))
clusterEvalQ(c1, library(boot))
clusterEvalQ(c1, library(data.table))
clusterEvalQ(c1, library(Rfast))
clusterEvalQ(c1, library(dplyr))
clusterEvalQ(c1, library(rjags))
clusterEvalQ(c1, library(R2jags))
clusterEvalQ(c1, library(elrm))
clusterEvalQ(c1, library(BBmisc))

start = Sys.time()
BAR_results_JAGS_812      <-  parLapply(c1, 1:12, function(i){    # Evaluate 6 scenarios 
  replicate(niterations, 
            simulation_replicate_jags(alpha = sim_scenario[[i]]$alpha,  beta = sim_scenario[[i]]$beta, gamma = sim_scenario[[i]]$gamma,  delta= sim_scenario[[i]]$delta,
                                      biomarkerprevalences,K,a,b,maxN,n.eachstage,delay,recruitmentpermonth), simplify = FALSE)
})
stopCluster(c1)

stop = Sys.time()
print(paste("Run time is:", stop-start,sep=""))

# Save results as RDS file to a folder
saveRDS(BAR_results_JAGS_812, "Multiple Biomarkers project/results/BAR_results_JAGS_812")

BAR_results_JAGS_812 = BAR_results_JAGS

#--------------------------------------------------------------------------------------------------

BAR_JAGS_operationalxtics = c("teststatistics", "onesidedpvalues", "allocation", "response", "biomarkerprofiles", "allocationpersubgroup", "patientsontreatment", 
                              "patients_per_arm", "responders_per_arm", "patientsonBesttreatment")
operational_xtics_compiled <- lapply(1:length(BAR_results_JAGS_812), function(i) {
  
  lapply(1:length(BAR_JAGS_operationalxtics), function(j){
    
    sapply(BAR_results_JAGS_812[[i]], "[[", BAR_JAGS_operationalxtics[j])
  })
})

names(operational_xtics_compiled)        <- paste0(rep("scenario",2), c(1:12))

for (j in 1:length(operational_xtics_compiled)) {
  
  names(operational_xtics_compiled[[j]]) <- BAR_JAGS_operationalxtics 
}

teststatistics_sims <- lapply(operational_xtics_compiled, "[[", c('teststatistics'))

patients_per_arm_sims        <- lapply(operational_xtics_compiled, "[[", c('patients_per_arm'))
responders_per_arm_sims      <- lapply(operational_xtics_compiled, "[[", c('responders_per_arm'))
patientsontreatment_sims     <- lapply(operational_xtics_compiled, "[[", c('patientsontreatment')) 
patientsonBestTreatment_sims <- lapply(BAR_results_JAGS_812, "[[", c('patientsonBestTreatment')) # Haven't extracte dthis yet in the code

# Statistical power
powerteststatistics <- sapply(1:length(teststatistics_sims),function(x){apply(teststatistics_sims[[x]], 1, function(i){sum(i> 1.645)/niterations})})
rownames(powerteststatistics) <- c("treatment1", "treatment1:B1", "treatment2:B1","treatment3:B1", "treatment4:B1",	
                                   "treatment2", "treatment1:B2", "treatment2:B2", "treatment3:B2","treatment4:B2",
                                   "treatment3", "treatment1:B3", "treatment2:B3", "treatment3:B3", "treatment4:B3", 
                                   "treatment4",	"treatment1:B4", "treatment2:B4",	"treatment3:B4"	,"treatment4:B4")

# 3.  On average how many patients go to each arm
summary_treatment_allocation <- sapply(1:length(patients_per_arm_sims),function(x){apply(patients_per_arm_sims[[x]], 1, summary)})

# 4. On average how any patients respond in each Tx arm
summary_treatment_nonresponse <- sapply(1:length(responders_per_arm_sims),function(x){ apply(responders_per_arm_sims[[x]][c(T, F),], 1, summary)}) # non-responders/0 Outcome
summary_treatment_response    <- sapply(1:length(responders_per_arm_sims),function(x){ apply(responders_per_arm_sims[[x]][c(F, T),], 1, summary)})# responders/1 Outcome

# 5. On average how many patients are on treatment
summary_patientsontreatment   <- sapply(1:length(patientsontreatment_sims),function(x){apply(rbind(patientsontreatment_sims[[x]]), 1, summary)})

# # 6. How many patients (receive)allocated to the best treatment
# #   Proportion of patients that get the best Tx available for them
summary_patientsonBestTreatment   <- sapply(1:length(patientsonBestTreatment_sims),function(x){apply(rbind(patientsonBestTreatment_sims[[x]]), 1, summary)})

#------------------ Export simulation results to excel worksheet ----------------------------------------------------------#
operational_xtics_combined <- list(powerteststatistics=powerteststatistics, summary_treatment_allocation=summary_treatment_allocation, 
                                   summary_treatment_response=summary_treatment_response, summary_treatment_nonresponse=summary_treatment_nonresponse,
                                   summary_patientsontreatment=summary_patientsontreatment, summary_patientsonBestTreatment=summary_patientsonBestTreatment)

BAR_JAGS_operationalxtics_final   <- names(operational_xtics_combined)

for (j in 1:length(operational_xtics_combined)) {
  
  colnames(operational_xtics_combined[[j]]) <- paste0(rep("scenario",2), c(1:12))
}     

# Export results to differnet worksheets of the same file
wb <- createWorkbook()

for (j in 1:length(BAR_JAGS_operationalxtics_final)) {
  addWorksheet(wb, BAR_JAGS_operationalxtics_final[j])
  writeData(wb=wb, sheet = BAR_JAGS_operationalxtics_final[j], x = t(operational_xtics_combined[[j]]), rowNames = T, colNames = T, borders = "columns")
}

saveWorkbook(wb, file = "BAR_results_JAGS_812.xlsx", overwrite = T)
# ------------------------------------------------------------------------------
#                         END OF CODE 
#-------------------------------------------------------------------------------