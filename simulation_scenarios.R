##### File #########################################################################
# Author: Luke Ouma
# Last edited: 07/05/2020
# Details: File to define the simulation scenarios


# alpha-effect of control Tx
alpha1 = -1.1

# beta-effect of treatment 
beta1 = c(0,0,0,0)
beta2 = c(0.25,0,0,0)

# gamma- effect of biomarker
gamma1 = rep(0,4)

# delta- treatment-biomarker interaction effect
delta1 = matrix(0,4,4)
delta2 = delta1; delta2[1,1]= 1.32 
delta3 = delta1; delta3[1,2]= 1.32  # T1 provides benefit for B2+ patients
delta4 = delta1; delta4[1,1]= -1.32  # T1 has detrimental effect for B1+ patients
delta5 = delta1; delta5[1,2]= -1.32  # T1 has detrimental effect for B2+ patients
delta6 = delta1; delta6[1,1]= 1.32; delta6[1,2]= 1.32 # T1 benefits 2 subgroups B1+ and B2+
delta7 = delta1; delta7[1,1]= -1.32; delta7[1,2]= -1.32 # T1 harms 2 subgroups B1+ and B2+
delta8 = delta1; delta8[1,1]= 1.32; delta8[1,2]= -1.32 # T1 benefits B1+ and harm in B2+ 
delta9 = delta1; delta9[1,1]= -1.32; delta9[1,2]= 1.32 # T1 harms B1+ and benefit in B2+ 
delta10 = delta1; delta10[1,1]= -1.32; delta10[2,1]= 1.32 # T1 provides harm in B1+, T2 provides benefit in B1+ 


sim_scenario = list(scenario1  = list(alpha=alpha1,  beta=beta1,  gamma=gamma1,  delta=t(delta1)),
                    scenario2  = list(alpha=alpha1,  beta=beta1,  gamma=gamma1,  delta=t(delta2)),
                    scenario3  = list(alpha=alpha1,  beta=beta1,  gamma=gamma1,  delta=t(delta3)),
                    scenario4  = list(alpha=alpha1,  beta=beta1,  gamma=gamma1,  delta=t(delta4)) ,
                    scenario5  = list(alpha=alpha1,  beta=beta1,  gamma=gamma1,  delta=t(delta5)),
                    scenario6  = list(alpha=alpha1,  beta=beta1,  gamma=gamma1,  delta=t(delta6)),
                    scenario7  = list(alpha=alpha1,  beta=beta1,  gamma=gamma1,  delta=t(delta7)),
                    scenario8  = list(alpha=alpha1,  beta=beta1,  gamma=gamma1,  delta=t(delta8)),
                    scenario9  = list(alpha=alpha1,  beta=beta1,  gamma=gamma1,  delta=t(delta9)),
                    scenario10 = list(alpha=alpha1,  beta=beta1,  gamma=gamma1,  delta=t(delta10)),
                    
                    scenario11 = list(alpha=alpha1,  beta=beta2,  gamma=gamma1,  delta=t(delta1)),
                    scenario12 = list(alpha=alpha1,  beta=beta2,  gamma=gamma1,  delta=t(delta2)))

