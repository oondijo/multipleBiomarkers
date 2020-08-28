##### File #####################################################################
# Author: Luke Ouma
# Last edited: 29/04/2020
# Details: File to evaluate operational characteristics for randomisation, hierarchy
#          and constrained randomization (allocation to rarer group) approaches of Tx 
#          allocation in the presence of multiple biomarkers

################# load relevant packages ##########
if(!require("Rfast")) install.packages("Rfast"); library(Rfast); 
if(!require("data.table")) install.packages("data.table"); library(data.table)
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)

################ FUNCTION ######################################################
multiple_biomarker_allocation <- function(
  prevalence, npatients, 
  alpha,  beta, gamma,  delta,  allocation_rule, 
  prob_hierarchy, prob_control, prob_lower, I)   {
  
  set.seed(I*169254 + 7382)
  ##### Variable initialisation ################################################

  nbiomarkers                <- length(prevalence)
  ntreatments                <- nbiomarkers + 1L       # Control and 4 experimental treatments
  
  # Treatment eligiblity matrix; everyone is eligible for the control
  treatment_eligibility      <- matrix(0L, npatients, ntreatments)
  treatment_eligibility[, 1] <- 1L
  
  ##### Random biomarker generation
  
  # Biomarkers will be the covariates for which we want to balance. They solely determine which treatments a patient is 
  # eligible for.
  
  for (biomarker in 1:nbiomarkers) {
    treatment_eligibility[, 1 + biomarker] <-
      stats::rbinom(npatients, 1, prevalence[biomarker])
  }
  
  biomarkers          = treatment_eligibility[, -1]
  
  # Find those negative for all biomarkers and make them eligible for all treatments
  
  treatment_eligibility[which(rowSums(treatment_eligibility[, -1]) ==
                                0), ]      <- 1L
  
  
  ###### Treatment allocation
  
  # Call different functions for treatment allocation
  if (allocation_rule        == "randomization") {
    treatment_allocated <- equal_randomization_eligible_Tx(npatients = npatients,  
                                                           treatment_eligibility= treatment_eligibility)
  } else if (allocation_rule == "prioritize_randomization") {
    treatment_allocated <- prioritize_randomization_eligible_Tx (npatients = npatients, prob_control = prob_control, 
                                                                 treatment_eligibility = treatment_eligibility)
  } else if (allocation_rule == "hierarchy"){
    treatment_allocated <- hierarchy(npatients = npatients, prob_hierarchy = prob_hierarchy, 
                                     treatment_eligibility= treatment_eligibility)
  } else if (allocation_rule == "constrained_randomization") {
    treatment_allocated <- constrained_randomization(npatients = npatients, ntreatments = ntreatments, prob_lower  = prob_lower, 
                                                     treatment_eligibility= treatment_eligibility) 
  }
  
  
  ############## Simulate outcome ######################
  
  outcome_logit <- rep(0,npatients)
  outcome_logit <- alpha + c(0,beta)[treatment_allocated] + as.double(biomarkers%*%gamma)
  
  for(i in 1:npatients)
  {
    for(biomarker in 1:nbiomarkers)
    {
      outcome_logit[i] <- outcome_logit[i] + biomarkers[i,biomarker]*delta[treatment_allocated[i],biomarker]
      
    }
  }
  
  ###  Make Dummy function
  make_dummies  <- function(v, prefix = '') {
    s <- sort(unique(v))
    d <- outer(v, s, function(v, s) 1L * (v == s))
    colnames(d) <- paste0(prefix, s)
    d
  }
  
  # Probability of outcome given allocated treatment
  outcome_prob    <- exp(outcome_logit)/(1+ exp(outcome_logit))
  outcome         <- stats::rbinom(npatients, 1, prob = outcome_prob)
  
  data            <- cbind.data.frame(biomarkers, 
                                      make_dummies(treatment_allocated, prefix = "treatment"), 
                                      outcome_prob, outcome)
  setnames(data,c("1","2","3","4"),c("B1","B2", "B3", "B4"))
  
  # Number allocated to each arm
  patients_per_arm   <- table(treatment_allocated)
  
  # Number of responders on each arm
  responders_per_arm <- table(outcome, treatment_allocated)
  
  # Number/Proportion of patients being treated
  patientsontreatment <- (npatients - sum(data$treatment1)) / npatients *100
  
  # No patients on the best Treatment available to them
  #.....abit of data manipulation here.....
  data$treatment      <- numeric(length = npatients)
  dataTx              <- data %>% dplyr::select(starts_with("treatment")) 
  for (i in 1:nrow(dataTx)) {
    data$treatment[i] <- unname(which.max(apply(dataTx[i,], MARGIN=2, max)))
  } 
  
  biom_outcomeprob_treat <- data %>% group_by(B1, B2,B3,B4, outcome_prob, treatment) %>% summarise(count =n()) %>% tidyr::unite(biomarker, B1:B4, sep="")
  aggregateBiomarkergrp_outcomeprob <- aggregate (count ~ biomarker + outcome_prob, data = biom_outcomeprob_treat,  sum )
  dup                               <- duplicated(aggregateBiomarkergrp_outcomeprob$biomarker)
  aggregateBiomarkergrp_outcomeprob <- cbind(aggregateBiomarkergrp_outcomeprob, dup)
  
  aggregateBiomarkergrp_outcomeprob <- aggregateBiomarkergrp_outcomeprob %>% arrange(biomarker, -outcome_prob)
  aggregateBiomarkergrp_outcomeprob <- aggregateBiomarkergrp_outcomeprob %>% group_by(biomarker) %>% filter(outcome_prob == max(outcome_prob))
  patientsonBestTreatment           <- (sum(aggregateBiomarkergrp_outcomeprob$count[aggregateBiomarkergrp_outcomeprob$dup=="TRUE"]))/ npatients *100
  

  ####### Fit GLM model #########################################
  
  glm_model    <- try(glm(outcome ~    treatment2 +  treatment3 + treatment4  + treatment5 +  # treatment effects
                            B1 + B2 + B3 + B4 +                                               # Biomarker effects
                            B1*treatment2 + B1*treatment3 + B1*treatment4 + B1*treatment5 +   # Tx-marker interaction
                            B2*treatment2 + B2*treatment3 + B2*treatment4 + B2*treatment5 +
                            B3*treatment2 + B3*treatment3 + B3*treatment4 + B3*treatment5 +
                            B4*treatment2 + B4*treatment3 + B4*treatment4 + B4*treatment5,
                          family = "binomial", data = data), silent = T)
  
  glm_model_pvalues2sided <- coef(summary(glm_model))[,4]
  glm_model_coefficients <- coef(summary(glm_model))[,1]
  
  # Perform a one sided test for interactions
  treatment2B1 = c(0, 1, rep(0,7), 1, rep(0, 15))
  treatment3B1 = c(rep(0,2), 1, rep(0,7), 1, rep(0, 14))
  treatment4B1 = c(rep(0,3), 1, rep(0,7), 1, rep(0, 13))
  treatment5B1 = c(rep(0,4), 1, rep(0,7), 1, rep(0, 12))
  
  treatment2B2 = c(0, 1, rep(0,11), 1, rep(0, 11))
  treatment3B2 = c(rep(0,2), 1, rep(0,11), 1, rep(0, 10))
  treatment4B2 = c(rep(0,3), 1, rep(0,11), 1, rep(0, 9))
  treatment5B2 = c(rep(0,4), 1, rep(0,11), 1, rep(0, 8))
  
  treatment2B3 = c(0, 1, rep(0,15), 1, rep(0, 7))
  treatment3B3 = c(rep(0,2), 1, rep(0,15), 1, rep(0, 6))
  treatment4B3 = c(rep(0,3), 1, rep(0,15), 1, rep(0, 5))
  treatment5B3 = c(rep(0,4), 1, rep(0,15), 1, rep(0, 4))
  
  treatment2B4 = c(0, 1, rep(0,19), 1, rep(0, 3))
  treatment3B4 = c(rep(0,2), 1, rep(0,19), 1, rep(0, 2))
  treatment4B4 = c(rep(0,3), 1, rep(0,19), 1, 0)
  treatment5B4 = c(rep(0,4), 1, rep(0,19), 1)
  
  matrixinteractions = rbind(treatment2B1,treatment3B1, treatment4B1,treatment5B1,
                             treatment2B2, treatment3B2,treatment4B2 ,treatment5B2 ,treatment2B3 ,
                             treatment3B3,treatment4B3 ,treatment5B3 ,treatment2B4 ,treatment3B4,
                             treatment4B4 ,treatment5B4)
  
  
  # Check if any coefficients are NA
  whichna   <- is.na(as.double(glm_model$coef))
  varmatrix <- summary(glm_model)$cov.unscaled
  
  if(any(whichna)){
    
    #set corresponding coefficient to zero
    glm_model$coef[whichna] <- 0
    
    #Add a row and column to glm1$cov.unscaled with infinite entries
    tempmatrix <- matrix(0,25,25)
    tempmatrix[!whichna, !whichna] <- summary(glm_model)$cov.unscaled
    tempmatrix[whichna, ]          <- 9e50
    tempmatrix[, whichna]          <- 9e50
    varmatrix                      <- tempmatrix
  }
  
  # Check whether model converged
  convergence <- glm_model$converged
  
  pvaluesinteractions <- sapply(1:nrow(matrixinteractions), function(x) {
    pnorm(-abs((matrixinteractions[x,] %*% glm_model$coef) / sqrt(t(matrixinteractions[x,]) %*% varmatrix %*% matrixinteractions[x,])))
  })
  
  #Compute test statistic
  teststatistics               <- sapply(1:nrow(matrixinteractions), function(x) {
    (matrixinteractions[x,] %*% glm_model$coef) / sqrt(t(matrixinteractions[x,]) %*% varmatrix %*% matrixinteractions[x,])
  })
  names(pvaluesinteractions)  <- names(glm_model$coef[10:25]) 

  names(teststatistics)  <- names(glm_model$coef[10:25]) 
  
  # Keep elements to be returned in a list
  required_output= list(convergence= convergence, glm_model_pvalues2sided = glm_model_pvalues2sided, teststatistics=teststatistics,
                        pvaluesinteractions = pvaluesinteractions,patientsontreatment = patientsontreatment, 
                        patientsonBestTreatment = patientsonBestTreatment,
                        patients_per_arm = patients_per_arm, responders_per_arm = responders_per_arm)
  return(required_output)
}