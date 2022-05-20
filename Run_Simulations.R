library(dplyr) ;library(ggplot2) ; library(RxODE) ; library(pracma)

#Select appropriate working directory
#setwd("________")

#Runs script which samples parameters for simulation

source("./Simulate_Parameters.R")

#The models:

#PK model alone, for Cav calculation over first two administrations:

PK_access_target <- RxODE({
  k12 = Q/V1
  k21 = Q/V2
  ke = CL/V1
  d/dt(C) = -(ke + k12) * C + k21 * P/V1 - kon * C * sTYRP1 + koff * Complex
  d/dt(P) = +k12 * C*V1 - k21 * P
  
  d/dt(sTYRP1) = -kon * C * sTYRP1 + koff * Complex
  d/dt(Complex) = +kon * C * sTYRP1 - koff * Complex
  
  
  Cc_mgL = C * 194000/1e9 * 1000 #Change of units using molecular weight (nM -> mgL)
})


Delay_between_baseline_and_admin <- 7 #From baseline scan of patient until start to treatment, how much time in days?

#PD model:

TGI_model <- RxODE({
  #PK needed for killing Emax function
  k12 = Q/V1
  k21 = Q/V2
  ke = CL/V1
  d/dt(C) = -(ke + k12) * C + k21 * P/V1 - kon * C * sTYRP1 + koff * Complex
  d/dt(P) = +k12 * C*V1 - k21 * P
  
  d/dt(sTYRP1) = -kon * C * sTYRP1 + koff * Complex
  d/dt(Complex) = +kon * C * sTYRP1 - koff * Complex
  
  
  Cc_mgL = C * 194000/1e9 * 1000
  
  #Killing Emax function
  Kill_TCB = Kill_max * Cc_mgL/(IC50 + Cc_mgL)
  
  #Time_to_resistance as function of Cav_individual_i:
  
  Time_to_Regrowth = Max_Time_to_Regrowth * Cav^hill_TtR/(Cav^hill_TtR + Cav50_TtR^hill_TtR)
  
  #KL_resistant as function of Cav_individual_i:
  
  KL_resistant = Max_KL_resistant * Cav/(Cav + Cav50_KLres)
  
  # The TGI model
  if(time < Time_to_Effect){d/dt(TS) = KL_sensible * TS}
  else if(time > Time_to_Effect & time < Time_to_Regrowth){d/dt(TS) = (KL_sensible - Kill_TCB) * TS}
  else{d/dt(TS) = KL_resistant * TS}
  
})


#For combination, KLres and Time_to_Regrowth are explicitly defined, independnet
#from Caverage

TGI_model_combo <- RxODE({
  #PK needed for killing Emax function
  k12 = Q/V1
  k21 = Q/V2
  ke = CL/V1
  d/dt(C) = -(ke + k12) * C + k21 * P/V1 - kon * C * sTYRP1 + koff * Complex
  d/dt(P) = +k12 * C*V1 - k21 * P
  
  d/dt(sTYRP1) = -kon * C * sTYRP1 + koff * Complex
  d/dt(Complex) = +kon * C * sTYRP1 - koff * Complex
  
  
  Cc_mgL = C * 194000/1e9 * 1000
  
  #Killing Emax function
  Kill_TCB = Kill_max * Cc_mgL/(IC50 + Cc_mgL)
  
  # The TGI model
  if(time < Time_to_Effect){d/dt(TS) = KL_sensible * TS}
  else if(time > Time_to_Effect & time < Time_to_Regrowth){d/dt(TS) = (KL_sensible - Kill_TCB) * TS}
  else{d/dt(TS) = KL_resistant * TS}
  
})

#Run PK model to get Caverage over first two cycles to improve speed of simulations:
Dose_mg <- 150

#Run takes about 12 seconds

NN = nrow(Parameters)
Cav_150mg <- numeric(length = NN)
for(i in 1:NN){
  Params_i <- Parameters[i, ]
  Sim_times_PK <- eventTable()
  Sim_times_PK <- add.sampling(eventTable = Sim_times_PK, time = seq(0, 42, length.out = 1000))
  Sim_times_PK <- add.dosing(eventTable = Sim_times_PK, dose = (Dose_mg/194000 * 1e9/1000)/
                               as.numeric(Params_i[which(names(Params_i) == "V1")]), 
                             nbr.doses = 2, dosing.interval = 21, do.sampling = F, start.time = 0)
  
  Solve_PK <- rxSolve(PK_access_target, params = Params_i, 
                      inits = c(C = 0, P = 0, Complex = 0, 
                                sTYRP1 = as.numeric(Params_i[which(names(Params_i) == "sTYRP1")])),
                      events = Sim_times_PK)
  Cav_i <- trapz(x = Solve_PK$time, y = Solve_PK$Cc_mgL)/42
  Cav_150mg[i] <- Cav_i
}


Parameters$Cav <- Cav_150mg


###################################################################################
#############     Run simulations of TGI at 150 mg, monotherapy      ##############
###################################################################################

# KL = MEDIUM ; KILLING = PROPORTIONAL to GROWTH

#This takes about two minutes, could be improved with pre-allocation

Parameters$KL_sensible <- KL_sensible_medium ; Parameters$Kill_max <- Kill_max_prop_medium 
Parameters$Max_KL_resistant <- Max_KL_resistant_medium


Results_150mg_medium_KL_prop_killing <- data.frame()
for(i in 1:NN){
  Params_i <- Parameters[i, ]
  Sim_times_PD <- eventTable()
  Sim_times_PD <- add.sampling(eventTable = Sim_times_PD, time = seq(0, 400, length.out = 1000))
  Sim_times_PD <- add.dosing(eventTable = Sim_times_PD, dose = (Dose_mg/194000 * 1e9/1000)/
                               as.numeric(Params_i[which(names(Params_i) == "V1")]), nbr.doses = 19, 
                             dosing.interval = 21, do.sampling = T, start.time = 0 + Delay_between_baseline_and_admin)
  
  Solve_PD <- rxSolve(TGI_model, params = Params_i, 
                      inits = c(C = 0, P = 0, Complex = 0, 
                                sTYRP1 = as.numeric(Params_i[which(names(Params_i) == "sTYRP1")]), TS = 50),
                      events = Sim_times_PD)
  Solve_PD$Individual_NN <- i
  Results_150mg_medium_KL_prop_killing <- rbind(Results_150mg_medium_KL_prop_killing, Solve_PD)
}



###################################################################################
#############     Run simulations of TGI at 150 mg, combination      ##############
###################################################################################

# KL = MEDIUM ; KILLING = PROPORTIONAL
#This takes about two minutes, could be improved with pre-allocation

#Re-assign the parameters to the ones of combination group
Parameters$KL_sensible <- KL_sensible_medium ; Parameters$Kill_max <- Kill_max_prop_medium
Parameters$KL_resistant <- KL_resistant_medium_combo
Parameters$Time_to_Regrowth <- Time_to_Regrowth_combo


Results_150mg_medium_KL_prop_killing_combo <- data.frame()
for(i in 1:NN){
  Params_i <- Parameters[i, ]
  Sim_times_PD <- eventTable()
  Sim_times_PD <- add.sampling(eventTable = Sim_times_PD, time = seq(0, 400, length.out = 1000))
  Sim_times_PD <- add.dosing(eventTable = Sim_times_PD, dose = (Dose_mg/194000 * 1e9/1000)/as.numeric(Params_i[which(names(Params_i) == "V1")]), nbr.doses = 19, dosing.interval = 21, do.sampling = T, start.time = 0 + Delay_between_baseline_and_admin)
  
  Solve_PD <- rxSolve(TGI_model_combo, params = Params_i, inits = c(C = 0, P = 0, Complex = 0, sTYRP1 = as.numeric(Params_i[which(names(Params_i) == "sTYRP1")]), TS = 50), events = Sim_times_PD)
  Solve_PD$Individual_NN <- i
  Results_150mg_medium_KL_prop_killing_combo <- rbind(Results_150mg_medium_KL_prop_killing_combo, Solve_PD)
}

