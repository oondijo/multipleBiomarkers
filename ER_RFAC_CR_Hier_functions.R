##### File #########################################################################
# Author: Luke Ouma
# Last edited: 
# Details: File to define the treatment allocation functions: 
#         -i)   Equal Randomization
#         -ii)  Randomization with fixed allocation probability to control
#         -iii) Hierarchy of biomarkers 
#         -iv)  Constrained randomization (allocation to the rarer group)
# 

###################################################################################
#         1.    Randomization 
#
################## Function #######################################################
equal_randomization_eligible_Tx <- function(
  npatients,                                   # No of patients
  treatment_eligibility) {    
  
  ##### Assign treatments ######################################################
  
  treatment_allocated           <- integer(npatients)
  
  for(patient in 1:npatients)
  {
    prob                        <- c(treatment_eligibility[patient,]/sum(treatment_eligibility[patient,])) 
    
    treatment_allocated[patient] <- sample(1:5,prob = prob,replace=T,size=1)
    
  }
  
  return(treatment_allocated)
}

#########################################################################################################
#  2. Randomization with fixed probability of allocation to control

################## Function #######################################################
prioritize_randomization_eligible_Tx <- function(
  npatients,                                  
  prob_control, treatment_eligibility)  {               
  
  ##### Assign treatments #######################################################
  
  treatment_allocated            <- integer(npatients)
  
  for(patient in 1:npatients)
  {

    prob                         <- c(prob_control, 
                                      (1 - prob_control)*treatment_eligibility[patient,-1]/sum(treatment_eligibility[patient,-1]))
    
    treatment_allocated[patient] <- sample(1:5, prob = prob, replace=T, size=1)
    
  }
  
  return(treatment_allocated)
}

#############################################################################################################
#              3.  Hierarchy of biomarkers
#############################################################################################################
hierarchy <- function(
  npatients,         
  prob_hierarchy,    
  treatment_eligibility)                 
  
{ 
  ##### Assign treatments #################################################
  
  treatment_allocated      <- integer(npatients)
  for (j in 1:npatients) {
    
    # Which treatments can a patient be allocated to
    treat                  <- which(treatment_eligibility[j,]==1)
    ntreat                 <- length(treat)
    
    next_treat             <- seq(treat)[1:2] # 1 is control, 2 is first treatment of the hierarchy a pt is eligible for. 
                                              # Treatment order always ascending. For hierarchy, assume 2 is the most important
    
    len_next_treat         <- length(next_treat)
    
    #### Allocation probabilities 
    # i) If a pt is eligible for 2 treatments only, i.e control and 1 marker-linked treatment,randomise between these two.
    # ii) If eligible for >1 treatment, allocate to either control or most important marker.
    #     Here, can define a probability of allocation to that important marker
    
    if (ntreat == 2) {   # i.e eligible for control and 1 marker-linked treatment only
      prob                 <- rep(1, ntreat)/ntreat
    } else {prob           <-  numeric(ntreat)
    
    # Randomise equally btn control and most important marker of the hierarchy in the case of multiple biomarkers
    
    prob[next_treat]       <- prob_hierarchy/len_next_treat
    prob[-next_treat]      <- (1 - prob_hierarchy)/(ntreat - len_next_treat)
    }
    treatment_allocated[j] <- sample(treat, prob = prob, replace = T, size = 1)
    
  }  
  
  return(treatment_allocated)
}

##############################################################################################################
#                 4. Constrained Randomization (Allocation to the rarer subgroup)
##############################################################################################################

constrained_randomization  <- function(
  npatients, ntreatments,
  prob_lower,                    # Probability allocated to Tx arms with lower patient accrual
  treatment_eligibility) 
  
{      
  # Assign treatments
  
  treatment_allocated      <- integer(npatients)
  treatment_allocated[1]   <- sample(1:ntreatments, 1, replace = T,
                                     prob = treatment_eligibility[1, ]/
                                       sum(treatment_eligibility[1, ]))
  for (j in 2:npatients) {
    
    # check Which treatments is this patient eligible for
    treat                  <- which(treatment_eligibility[j, ] == 1)
    ntreat                 <- length(treat)
    
    # Check how many patients already allocated to treatments for which next patient is eligible
    patients_treat         <- integer(ntreat)
    for (i in 1:ntreat) {
      patients_treat[i]    <- sum(treatment_allocated[1:(j - 1)] == treat[i])
    }
    
    # Check which Treatments have the lowest allocation so far
    next_treat             <- which(patients_treat == min(patients_treat))
    
    # Determine next treatment allocation probabilities. Two cases can occur:
    # i)  either equal no of pts already allocated to ALL treatments -> equal chance of
    #     allocation to all
    # ii) at least one of the treatments has more/fewer patients -> 100*prob_lower% chance
    #     allocation to treatment with fewer pts
    
    len_nexttreat          <- length(next_treat)
    if (len_nexttreat == ntreat) {
      prob                 <- rep(1, ntreat)/ntreat
    } else {
      prob                 <- numeric(ntreat)
      prob[next_treat]     <- prob_lower/len_nexttreat
      prob[-next_treat]    <- (1 - prob_lower)/(ntreatments - len_nexttreat)
    }
    
    # Allocate treatment
    treatment_allocated[j] <- sample(treat, prob = prob, replace = T, size = 1)
  }
  
  return(treatment_allocated)
  
}

###################################################################################################################