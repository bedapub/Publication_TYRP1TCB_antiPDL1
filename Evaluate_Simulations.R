library(ggplot2) ; library(pracma) ; library(dplyr)

source("./Run_Simulations.R")

######################################################################################
###########     Evaluations of responses with simulations and scans       ############
######################################################################################


#Use the full granularity of the simulation to derive clinical metrics of response

#Idea is as follows: calculate nadir (min TS value during simulation), check if
# at some point we have at least 30% reduction in SLD (ie, TS<35) -> response,
# then check when we have a 20% increase from nadir -> progression
# Time difference between first response and first
# progression is DoR; time until first progression is PFS


Evaluate_150mg_medium_KL_prop_killing <- Results_150mg_medium_KL_prop_killing %>% group_by(Individual_NN) %>%
  summarise(NADIR = min(TS), NADIR_time = first(time[which(TS == min(TS))]),
            Start_Response = if(sum(TS<35)>=1){first(time[which(TS < 35)])} else{NA},
            Start_Progression = if(is.na(first(time[which(TS > 1.2 * min(TS) & time > time[first(which(TS == min(TS)))])]))){max(Results_150mg_medium_KL_prop_killing$time)} else{first(time[which(TS > 1.2 * min(TS) & time > time[first(which(TS == min(TS)))])])},
            Time_below_baseline_progression = if(is.na(first(time[which(TS > 1.2 * 50)]))){max(Results_150mg_medium_KL_prop_killing$time)} else{first(time[which(TS > 1.2 * 50)])},
            Duration_of_Response = if(is.na(first(time[which(TS > 1.2 * min(TS) & time > time[first(which(TS == min(TS)))])]))){max(Results_150mg_medium_KL_prop_killing$time)} 
            else if(is.na(Start_Response)){NA} else{Start_Progression  - Start_Response})


Evaluate_150mg_medium_KL_prop_killing_combo <- Results_150mg_medium_KL_prop_killing_combo %>% group_by(Individual_NN) %>%
  summarise(NADIR = min(TS), NADIR_time = first(time[which(TS == min(TS))]),
            Start_Response = if(sum(TS<35)>=1){first(time[which(TS < 35)])} else{NA},
            Start_Progression = if(is.na(first(time[which(TS > 1.2 * min(TS) & time > time[first(which(TS == min(TS)))])]))){max(Results_150mg_medium_KL_prop_killing_combo$time)} else{first(time[which(TS > 1.2 * min(TS) & time > time[first(which(TS == min(TS)))])])},
            Time_below_baseline_progression = if(is.na(first(time[which(TS > 1.2 * 50)]))){max(Results_150mg_medium_KL_prop_killing_combo$time)} else{first(time[which(TS > 1.2 * 50)])},
            Duration_of_Response = if(is.na(first(time[which(TS > 1.2 * min(TS) & time > time[first(which(TS == min(TS)))])]))){max(Results_150mg_medium_KL_prop_killing_combo$time)} 
            else if(is.na(Start_Response)){NA} else{Start_Progression  - Start_Response})


##########################################################################
###########  RECIST 1.1 evaluation with scans every 6 weeks   ############
##########################################################################

#We select TS measurements every six weeks
Desired_times_six <- c(0, seq(from = 6*7 + Delay_between_baseline_and_admin, to = 400, by = 6*7))

#In our simulation, which are the closest time points to the ones we want=
Scan_times_six <- Results_150mg_medium_KL_prop_killing$time[
  sapply(Desired_times_six, function(x) 
    which.min(abs(Results_150mg_medium_KL_prop_killing$time - x)))]

#Filter the TS at the desired times
Scans_6_weeks_150mg_medium_KL_prop_killing <- Results_150mg_medium_KL_prop_killing %>% 
  filter(time %in% Scan_times_six)
Scans_6_weeks_150mg_medium_KL_prop_killing_combo <- Results_150mg_medium_KL_prop_killing_combo %>% 
  filter(time %in% Scan_times_six)


#Evaluate with RECIST: 

Evaluate_Scans_6_weeks_150mg_medium_KL_prop_killing <-	Scans_6_weeks_150mg_medium_KL_prop_killing %>% group_by(Individual_NN) %>%
  summarise(NADIR = min(TS), NADIR_time = first(time[which(TS == NADIR)]),
            Start_Response = if(sum(TS<35)>=1){first(time[which(TS < 35)])} else{NA},
            Start_Progression = if(is.na(first(time[which(TS > 1.2 * NADIR &            #Criteria for progression
                                                          time > time[first(which(TS == NADIR))])]))){ #No progression? Then max(sim time)
              max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{ #Progression? First time point which fullfils criteria
                first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])])},
            Time_below_baseline_progression = if(is.na(first(time[which(TS > 1.2 * 50)]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{first(time[which(TS > 1.2 * 50)])},
            Duration_of_Response = if(is.na(first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} 
            else if(is.na(Start_Response)){NA} else{Start_Progression  - Start_Response}) 

Evaluate_Scans_6_weeks_150mg_medium_KL_prop_killing_combo <-	Scans_6_weeks_150mg_medium_KL_prop_killing_combo %>% group_by(Individual_NN) %>%
  summarise(NADIR = min(TS), NADIR_time = first(time[which(TS == NADIR)]),
            Start_Response = if(sum(TS<35)>=1){first(time[which(TS < 35)])} else{NA},
            Start_Progression = if(is.na(first(time[which(TS > 1.2 * NADIR &            #Criteria for progression
                                                          time > time[first(which(TS == NADIR))])]))){ #No progression? Then max(sim time)
              max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{ #Progression? First time point which fullfils criteria
                first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])])},
            Time_below_baseline_progression = if(is.na(first(time[which(TS > 1.2 * 50)]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{first(time[which(TS > 1.2 * 50)])},
            Duration_of_Response = if(is.na(first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} 
            else if(is.na(Start_Response)){NA} else{Start_Progression  - Start_Response}) 

##########################################################################
###########  RECIST 1.1 evaluation with scans every 12 weeks   ############
##########################################################################

Desired_times_twelve <- c(0, seq(from = 12*7 + Delay_between_baseline_and_admin, to = 400, by = 12*7))
Scan_times_twelve <- Results_150mg_medium_KL_prop_killing$time[sapply(Desired_times_twelve, function(x) which.min(abs(Results_150mg_medium_KL_prop_killing$time - x)))]
Scans_12_weeks_150mg_medium_KL_prop_killing <- Results_150mg_medium_KL_prop_killing %>% filter(time %in% Scan_times_twelve)
Scans_12_weeks_150mg_medium_KL_prop_killing_combo <- Results_150mg_medium_KL_prop_killing_combo %>% filter(time %in% Scan_times_twelve)

Evaluate_Scans_12_weeks_150mg_medium_KL_prop_killing <-	Scans_12_weeks_150mg_medium_KL_prop_killing %>% group_by(Individual_NN) %>%
  summarise(NADIR = min(TS), NADIR_time = first(time[which(TS == NADIR)]),
            Start_Response = if(sum(TS<35)>=1){first(time[which(TS < 35)])} else{NA},
            Start_Progression = if(is.na(first(time[which(TS > 1.2 * NADIR &            #Criteria for progression
                                                          time > time[first(which(TS == NADIR))])]))){ #No progression? Then max(sim time)
              max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{ #Progression? First time point which fullfils criteria
                first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])])},
            Time_below_baseline_progression = if(is.na(first(time[which(TS > 1.2 * 50)]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{first(time[which(TS > 1.2 * 50)])},
            Duration_of_Response = if(is.na(first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} 
            else if(is.na(Start_Response)){NA} else{Start_Progression  - Start_Response}) 

Evaluate_Scans_12_weeks_150mg_medium_KL_prop_killing_combo <-	Scans_12_weeks_150mg_medium_KL_prop_killing_combo %>% group_by(Individual_NN) %>%
  summarise(NADIR = min(TS), NADIR_time = first(time[which(TS == NADIR)]),
            Start_Response = if(sum(TS<35)>=1){first(time[which(TS < 35)])} else{NA},
            Start_Progression = if(is.na(first(time[which(TS > 1.2 * NADIR &            #Criteria for progression
                                                          time > time[first(which(TS == NADIR))])]))){ #No progression? Then max(sim time)
              max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{ #Progression? First time point which fullfils criteria
                first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])])},
            Time_below_baseline_progression = if(is.na(first(time[which(TS > 1.2 * 50)]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{first(time[which(TS > 1.2 * 50)])},
            Duration_of_Response = if(is.na(first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} 
            else if(is.na(Start_Response)){NA} else{Start_Progression  - Start_Response}) 

##########################################################################
###########  RECIST 1.1 evaluation with scans every 9 weeks   ############
##########################################################################

Desired_times_nine <- c(0, seq(from = 9*7 + Delay_between_baseline_and_admin, to = 400, by = 9*7)) 
Scan_times_nine <- Results_150mg_medium_KL_prop_killing$time[sapply(Desired_times_nine, function(x) which.min(abs(Results_150mg_medium_KL_prop_killing$time - x)))]
Scans_9_weeks_150mg_medium_KL_prop_killing <- Results_150mg_medium_KL_prop_killing %>% filter(time %in% Scan_times_nine)
Scans_9_weeks_150mg_medium_KL_prop_killing_combo <- Results_150mg_medium_KL_prop_killing_combo %>% filter(time %in% Scan_times_nine)

Evaluate_Scans_9_weeks_150mg_medium_KL_prop_killing <-	Scans_9_weeks_150mg_medium_KL_prop_killing %>% group_by(Individual_NN) %>%
  summarise(NADIR = min(TS), NADIR_time = first(time[which(TS == NADIR)]),
            Start_Response = if(sum(TS<35)>=1){first(time[which(TS < 35)])} else{NA},
            Start_Progression = if(is.na(first(time[which(TS > 1.2 * NADIR &            #Criteria for progression
                                                          time > time[first(which(TS == NADIR))])]))){ #No progression? Then max(sim time)
              max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{ #Progression? First time point which fullfils criteria
                first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])])},
            Time_below_baseline_progression = if(is.na(first(time[which(TS > 1.2 * 50)]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{first(time[which(TS > 1.2 * 50)])},
            Duration_of_Response = if(is.na(first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} 
            else if(is.na(Start_Response)){NA} else{Start_Progression  - Start_Response}) 

Evaluate_Scans_9_weeks_150mg_medium_KL_prop_killing_combo <-	Scans_9_weeks_150mg_medium_KL_prop_killing_combo %>% group_by(Individual_NN) %>%
  summarise(NADIR = min(TS), NADIR_time = first(time[which(TS == NADIR)]),
            Start_Response = if(sum(TS<35)>=1){first(time[which(TS < 35)])} else{NA},
            Start_Progression = if(is.na(first(time[which(TS > 1.2 * NADIR &            #Criteria for progression
                                                          time > time[first(which(TS == NADIR))])]))){ #No progression? Then max(sim time)
              max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{ #Progression? First time point which fullfils criteria
                first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])])},
            Time_below_baseline_progression = if(is.na(first(time[which(TS > 1.2 * 50)]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} else{first(time[which(TS > 1.2 * 50)])},
            Duration_of_Response = if(is.na(first(time[which(TS > 1.2 * NADIR & time > time[first(which(TS == NADIR))])]))){max(Scans_6_weeks_150mg_medium_KL_prop_killing$time)} 
            else if(is.na(Start_Response)){NA} else{Start_Progression  - Start_Response}) 



#################################################################################################
###########     Best Overall Response, monotherapy and combination, with simulations ############
###########                        and different scanning patterns                   ############
#################################################################################################

# Values will be given in percentages, we divide by 10 because we have /1000 virtual patients * 100 to percentage = /10
BOR_rate_150mg_medium_KL_prop_killing <- sum(!(is.na(Evaluate_150mg_medium_KL_prop_killing$Start_Response)))/10


BOR_rate_150mg_medium_KL_prop_killing_combo <- sum(!(is.na(Evaluate_150mg_medium_KL_prop_killing_combo$Start_Response)))/10


#Requiring confirmation of the response in a subsequent scan:

Confirm_PR_Scans_6_weeks_150mg_medium_KL_prop_killing <- Scans_6_weeks_150mg_medium_KL_prop_killing %>% group_by(Individual_NN) %>%
  summarise(Scans_under_response_threshold = length(which(TS<35)))

Confirm_PR_Scans_6_weeks_150mg_medium_KL_prop_killing_combo <- Scans_6_weeks_150mg_medium_KL_prop_killing_combo %>% group_by(Individual_NN) %>%
  summarise(Scans_under_response_threshold = length(which(TS<35)))

Confirm_PR_Scans_12_weeks_150mg_medium_KL_prop_killing <- Scans_12_weeks_150mg_medium_KL_prop_killing %>% group_by(Individual_NN) %>%
  summarise(Scans_under_response_threshold = length(which(TS<35)))

Confirm_PR_Scans_12_weeks_150mg_medium_KL_prop_killing_combo <- Scans_12_weeks_150mg_medium_KL_prop_killing_combo %>% group_by(Individual_NN) %>%
  summarise(Scans_under_response_threshold = length(which(TS<35)))

Confirm_PR_Scans_9_weeks_150mg_medium_KL_prop_killing <- Scans_9_weeks_150mg_medium_KL_prop_killing %>% group_by(Individual_NN) %>%
  summarise(Scans_under_response_threshold = length(which(TS<35)))

Confirm_PR_Scans_9_weeks_150mg_medium_KL_prop_killing_combo <- Scans_9_weeks_150mg_medium_KL_prop_killing_combo %>% group_by(Individual_NN) %>%
  summarise(Scans_under_response_threshold = length(which(TS<35)))


ORR_6_mono <- length(which(Confirm_PR_Scans_6_weeks_150mg_medium_KL_prop_killing$
                             Scans_under_response_threshold >1))/10


ORR_6_combo <- length(which(Confirm_PR_Scans_6_weeks_150mg_medium_KL_prop_killing_combo$
                              Scans_under_response_threshold > 1))/10


ORR_12_mono <- length(which(Confirm_PR_Scans_12_weeks_150mg_medium_KL_prop_killing$
                              Scans_under_response_threshold > 1))/10


ORR_12_combo <- length(which(Confirm_PR_Scans_12_weeks_150mg_medium_KL_prop_killing_combo$
                               Scans_under_response_threshold > 1))/10


ORR_9_mono <- length(which(Confirm_PR_Scans_9_weeks_150mg_medium_KL_prop_killing$
                             Scans_under_response_threshold > 1))/10


ORR_9_combo <- length(which(Confirm_PR_Scans_9_weeks_150mg_medium_KL_prop_killing_combo$
                              Scans_under_response_threshold > 1))/10

######  Duration of response  #####

#Monotherapy

DoR_Sims_mono <- quantile(Evaluate_150mg_medium_KL_prop_killing$Duration_of_Response,
                          probs = c(0.05, 0.5, 0.95), na.rm = T)

DoR_6_mono <- quantile(Evaluate_Scans_6_weeks_150mg_medium_KL_prop_killing$Duration_of_Response,
                       probs = c(0.05, 0.5, 0.95), na.rm = T)

DoR_9_mono <- quantile(Evaluate_Scans_9_weeks_150mg_medium_KL_prop_killing$Duration_of_Response,
                       probs = c(0.05, 0.5, 0.95), na.rm = T)

DoR_12_mono <- quantile(Evaluate_Scans_12_weeks_150mg_medium_KL_prop_killing$Duration_of_Response,
                        probs = c(0.05, 0.5, 0.95), na.rm = T)


#Combination
DoR_Sims_combo <- quantile(Evaluate_150mg_medium_KL_prop_killing_combo$Duration_of_Response,
                           probs = c(0.05, 0.5, 0.95), na.rm = T)

DoR_6_combo <- quantile(Evaluate_Scans_6_weeks_150mg_medium_KL_prop_killing_combo$Duration_of_Response,
                        probs = c(0.05, 0.5, 0.95), na.rm = T)

DoR_9_combo <- quantile(Evaluate_Scans_9_weeks_150mg_medium_KL_prop_killing_combo$Duration_of_Response,
                        probs = c(0.05, 0.5, 0.95), na.rm = T)

DoR_12_combo <- quantile(Evaluate_Scans_12_weeks_150mg_medium_KL_prop_killing_combo$Duration_of_Response,
                         probs = c(0.05, 0.5, 0.95), na.rm = T)


######  PFS  #####

#Monotherapy

PFS_Sims_mono <- quantile(Evaluate_150mg_medium_KL_prop_killing$Start_Progression,
                          probs = c(0.05, 0.5, 0.95), na.rm = T)

PFS_6_mono <- quantile(Evaluate_Scans_6_weeks_150mg_medium_KL_prop_killing$Start_Progression,
                       probs = c(0.05, 0.5, 0.95), na.rm = T)

PFS_9_mono <- quantile(Evaluate_Scans_9_weeks_150mg_medium_KL_prop_killing$Start_Progression,
                       probs = c(0.05, 0.5, 0.95), na.rm = T)

PFS_12_mono <- quantile(Evaluate_Scans_12_weeks_150mg_medium_KL_prop_killing$Start_Progression,
                        probs = c(0.05, 0.5, 0.95), na.rm = T)


#Combination
PFS_Sims_combo <- quantile(Evaluate_150mg_medium_KL_prop_killing_combo$Start_Progression,
                           probs = c(0.05, 0.5, 0.95), na.rm = T)

PFS_6_combo <- quantile(Evaluate_Scans_6_weeks_150mg_medium_KL_prop_killing_combo$Start_Progression,
                        probs = c(0.05, 0.5, 0.95), na.rm = T)

PFS_9_combo <- quantile(Evaluate_Scans_9_weeks_150mg_medium_KL_prop_killing_combo$Start_Progression,
                        probs = c(0.05, 0.5, 0.95), na.rm = T)

PFS_12_combo <- quantile(Evaluate_Scans_12_weeks_150mg_medium_KL_prop_killing_combo$Start_Progression,
                         probs = c(0.05, 0.5, 0.95), na.rm = T)



Sim_mono = c(BOR_rate_150mg_medium_KL_prop_killing, DoR_Sims_mono, PFS_Sims_mono)
Scans_6_mono = c(ORR_6_mono, DoR_6_mono, PFS_6_mono)
Scans_9_mono = c(ORR_9_mono, DoR_9_mono, PFS_9_mono)
Scans_12_mono = c(ORR_12_mono, DoR_12_mono, PFS_12_mono)
Sim_combo = c(BOR_rate_150mg_medium_KL_prop_killing_combo, DoR_Sims_combo, PFS_Sims_combo)
Scans_6_combo = c(ORR_6_combo, DoR_6_combo, PFS_6_combo)
Scans_9_combo = c(ORR_9_combo, DoR_9_combo, PFS_9_combo)
Scans_12_combo = c(ORR_12_combo, DoR_12_combo, PFS_12_combo)

Response_Clinics <- as.data.frame(rbind(Sim_mono, Scans_6_mono, Scans_9_mono, Scans_12_mono,
                                        Sim_combo, Scans_6_combo, Scans_9_combo, Scans_12_combo))

names(Response_Clinics) <- c("ORR", "DoR_5th", "DoR_median", "DoR_95th", 
                             "PFS_5th", "PFS_median", "PFS_95th")


