##### File #########################################################################
# Author: Luke  Ouma
# Last edited: 29/04/2020
# Details: File to define the multiple biomarker superfunction
#         -i)   Run the simulations using all the approaches
#         -ii)  Compile and export all the results to workbooks (4 different workbooks, 
#               each with xx no of worksheets (xx is no of operational xtics interested in))
#   
###################################################################################

start = Sys.time()

##### Load scripts for the treatment allocation functions and simulation scenarios 
source("C:/Users/Luke Ondijo/OneDrive - Newcastle University/PhD research work/Multiple Biomarkers project/Rand_hierarchy_constrRand_functions.R")
source("C:/Users/Luke Ondijo/OneDrive - Newcastle University/PhD research work/Multiple Biomarkers project/operational_characteristics_rand_hier_constRand.R")
source("C:/Users/Luke Ondijo/OneDrive - Newcastle University/PhD research work/Multiple Biomarkers project/simulation_scenarios_multipleBiomarkers.R")

# Load relevant libraries
if(!require(foreach)) install.packages("foreach"); library(foreach)
if(!require(doParallel)) install.packages("doParallel"); library(doParallel)
if(!require(tidyverse)) install.packages("tidyverse"); library(tidyverse)
if(!require(openxlsx)) install.packages("openxlsx");library(openxlsx)

operating_xtics_super_fn <- function(
  prevalence, npatients, 
  alpha,  beta, gamma,  delta,  
  allocation_rule, reps,
  prob_hierarchy, prob_control, prob_lower, I) {
  
  cores  <- detectCores()
  c1     <- makeCluster(cores)
  registerDoParallel(c1)
  
  operating_xtics <- foreach(i = 1:reps, 
                             .export   = c('multiple_biomarker_allocation', 'constrained_randomization', 
                                           'equal_randomization_eligible_Tx', 'hierarchy', 'prioritize_randomization_eligible_Tx'),
                             .packages = c("boot","data.table","Rfast", "dplyr")) %dopar% {
                               
                               results = multiple_biomarker_allocation (
                                 alpha = alpha,  beta = beta, gamma = gamma,  delta= delta,
                                 allocation_rule = allocation_rule,
                                 prob_control = prob_control, prob_hierarchy = prob_hierarchy, prob_lower  = prob_lower,  
                                 npatients                  = npatients,         
                                 prevalence                 = prevalence,
                                 i)     
                               return(results)
                             }
  
  stopCluster(c1)
  
  ############# Operational xtics -summary ###########################################
  pvalues2sided_sims           <- lapply(operating_xtics, "[[", c('glm_model_pvalues2sided'))
  pvaluesinteractions_sims     <- lapply(operating_xtics, "[[", c('pvaluesinteractions'))
  patients_per_arm_sims        <- lapply(operating_xtics, "[[", c('patients_per_arm'))
  responders_per_arm_sims      <- lapply(operating_xtics, "[[", c('responders_per_arm'))
  patientsontreatment_sims     <- lapply(operating_xtics, "[[", c('patientsontreatment')) 
  patientsonBestTreatment_sims <- lapply(operating_xtics, "[[", c('patientsonBestTreatment')) 
  
  teststatistics_sims <- lapply(operating_xtics, "[[", c('teststatistics'))
  convergence_sims    <- lapply(operating_xtics, "[[", c('convergence')) 
  
  
  # 1. No of simulations where model converged # How many instances dont we have full model convergence
  sims_uncoverged           <- length(which(lengths(pvalues2sided_sims)<25))   
  prop_sims_uncoverged      <- sims_uncoverged/reps
  
  prop_modelsims_uncoverged <- length(which(lengths(convergence_sims)==FALSE))/reps   # True model convergence by checking
  
  # 2. Power - Check the probability that a treatment is better than control
   teststatistics_sims_2  <- do.call(rbind, teststatistics_sims)
  powerteststatistics    <- apply(teststatistics_sims_2, 2, function(x){sum(x> 1.645)/reps})
  
  
  # 3.  On average how many patients go to each arm
  patients_per_arm_sims2        <- do.call(rbind, patients_per_arm_sims)
  summary_treatment_allocation  <- apply(patients_per_arm_sims2, 2, summary)
  
  # 4. On average how any patients respond in each Tx arm
  responders_per_arm_sims2      <- do.call(rbind, responders_per_arm_sims)
  summary_treatment_nonresponse <- apply(responders_per_arm_sims2[c(T, F),], 2, summary)  # non-responders/0 Outcome
  summary_treatment_response    <- apply(responders_per_arm_sims2[c(F, T),], 2, summary)  # responders/1 Outcome
  
  # 5. On average how many patients are on treatment
  patientsontreatment_sims2     <- do.call(rbind, patientsontreatment_sims)
  summary_patientsontreatment   <- apply(patientsontreatment_sims2, 2, summary)
  
  # 6. How many patients (receive)allocated to the best treatment available for them
  patientsonBestTreatment_sims2     <- do.call(rbind, patientsonBestTreatment_sims)
  summary_patientsonBestTreatment   <- apply(patientsonBestTreatment_sims2, 2, summary)
  
  #################### Return results ##########
  operating_xtics_results = list(prop_sims_uncoverged,prop_modelsims_uncoverged,power_interaction,powerteststatistics,
                                 summary_treatment_allocation, 
                                 summary_treatment_nonresponse, summary_treatment_response,
                                 summary_patientsontreatment, summary_patientsonBestTreatment )
  return(operating_xtics_results)
}

allocation_rules = c("randomization", "prioritize_randomization", "hierarchy", "constrained_randomization")

mbiomarker_results = lapply(seq_along(allocation_rules), function(j){
  
  lapply(1:12, function(k){                  # Check the prespecified 12 scenarios
    operating_xtics_super_fn(
      npatients = 400, prevalence = c(0.3, 0.25, 0.3, 0.25),
      allocation_rule = allocation_rules[j], prob_control = 0.25, prob_hierarchy = 0.5 , prob_lower = 0.9,
      alpha = sim_scenario[[k]]$alpha,  beta = sim_scenario[[k]]$beta, gamma = sim_scenario[[k]]$gamma,  delta= sim_scenario[[k]]$delta,
      reps = 1e4)
  })
})

# Save results 
saveRDS(mbiomarker_results, "C:/Users/Luke Ondijo/OneDrive - Newcastle University/PhD research work/Multiple Biomarkers project/results/mbiomarker_results")

###########################################################################################################################

# Export simulation results to excel worksheet

setwd("C:/Users/Luke Ondijo/OneDrive - Newcastle University/PhD research work/Multiple Biomarkers project/results")

# load results
mbiomarker_results  <- readRDS("C:/Users/Luke Ondijo/OneDrive - Newcastle University/PhD research work/Multiple Biomarkers project/results/mbiomarker_results")

operational_xtics   <- c("prop_sims_uncoverged", "prop_modelsims_uncoverged","power_interaction", "powerteststatistics",
                         "summary_treatment_allocation", 
                         "summary_treatment_nonresponse", "summary_treatment_response", 
                         "summary_patientsontreatment", "summary_patientsonBestTreatment")


# Name/label each of different scenarios and results obtained from each scenario
for (j in 1:length(mbiomarker_results)) {
  
  names(mbiomarker_results[[j]])        <- paste0(rep("scenario",length(mbiomarker_results[[j]])), 1:length(mbiomarker_results[[j]]))
  
  for (k in 1:length(mbiomarker_results[[1]])) {   # there are identical no of sim scenarios for each Tx allocation procedure 
    
    names(mbiomarker_results[[j]][[k]]) <- operational_xtics
  } 
}

# Compile results for all 4 Treatment allocation approaches with all of the scenarios 
operational_xtics_compiled <- lapply(1:length(mbiomarker_results), function(i) {
  
  lapply(1:length(operational_xtics), function(j){
    
    sapply(mbiomarker_results[[i]], "[[", operational_xtics[j])
  })
})

# Export all results to 4 different excel files(Workbooks), each with 5 different worksheets
wb_list        <- list()
for (i in 1:length(allocation_rules)) {
  wb_list[[i]] <- createWorkbook()             # create a list of workbooks
  
  for (j in 1:length(operational_xtics)) {    # Create worksheets in each workbook for each x-tic we are interested in
    
    addWorksheet(wb_list[[i]], operational_xtics[j])
    
    writeData(wb = wb_list[[i]], sheet = operational_xtics[j], x = t(operational_xtics_compiled[[i]][[j]]), rowNames = T, colNames = T, borders="columns")
  }
  
  saveWorkbook(wb_list[[i]], file = paste(allocation_rules[i], ".xlsx", sep = ""), overwrite = T)
}


######################################### END #############################################################################     
stop = Sys.time()
print(paste("Run time is:", stop-start,sep=""))
