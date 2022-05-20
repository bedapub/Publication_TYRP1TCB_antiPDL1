# Simulate parameters from lognormal distr. to be used in TGI profile modeling

#Seeds for random number generation are changed for each parameter, to ensure
# no correlation at the individual parameter level is included


###   PK parameters escaled from mouse


NN = 1000 #Number of individuals to use in stochastic simulation

############
#    PK    #
#  params  #
############

# CL, clearance for TYRP1-TCB
set.seed(1)
Theta_CL = 1.093; Omega_CL = 0.275; CL = Theta_CL * exp(rnorm(n = NN, mean = 0, sd = Omega_CL))

# Q, intercompartmental clearance
set.seed(2)
Theta_Q = 3.826; Omega_Q = 0.1; Q = Theta_Q * exp(rnorm(n = NN, mean = 0, sd = Omega_Q))

#Volume of distribution of central compartment
set.seed(3)
Theta_V1 = 4.1279; Omega_V1 = 0.2; V1 = Theta_V1 * exp(rnorm(n = NN, mean = 0, sd = Omega_V1))

#Volume of distribution of peripheral compartment
set.seed(4)
Theta_V2 = 4.354; Omega_V2 = 0.1 ; V2 = Theta_V2 * exp(rnorm(n = NN, mean = 0, sd = Omega_V2)) 

#Accessible target from central compartment
set.seed(5)
sTYRP1 <- 0.9802707 * exp(rnorm(1000, 0, 1.465867))

#kon, koff, as function of Kd
Kd = 125 * 10^(-3) #Kd in nM
koff = 0.1 #assumed
kon = koff/Kd




############
##   PD   ##
## params ##
############

Delay_between_baseline_and_admin <- 7 #From baseline scan of patient until start to treatment, how much time in days?



#IC50, potency parameter of the killing function that relates TYRP1-TCB concentration to killing effect, with no variability
set.seed(6)
Theta_IC50 = 0.345; Omega_IC50 = 0.1 
IC50 = Theta_IC50 * exp(rnorm(n = NN, mean = 0, sd = Omega_IC50))


#Medium KL sensible from translational PKPD with pembrolizumab paper
set.seed(7)
Theta_KL_sensible_medium = 0.0036; Omega_KL_sensible = 0.52
KL_sensible_medium = Theta_KL_sensible_medium * exp(rnorm(n = NN, mean = 0, sd = Omega_KL_sensible))


#Kill_max, Emax of the killing function that relates TYRP1-TCB concentration to killing effect
#Proportional to growth constants (proportionality factor = Kill_max/KL_sens = 2.2975)

set.seed(10)
Theta_Kill_max_prop_medium = Theta_KL_sensible_medium * 2.2975; Omega_Kill_max = 0.27
Kill_max_prop_medium = Theta_Kill_max_prop_medium * exp(rnorm(n = NN, mean = 0, sd = Omega_Kill_max))


#Time to Effect, allometrically scaled
#We simulate 1000 values, then add to each of them the delay between baseline
#assessment and administration

set.seed(11)
Theta_Time_to_Effect_scaled = 9.7; Omega_Time_to_Effect = 0.64
Time_to_Effect_scaled = Theta_Time_to_Effect_scaled * exp(rnorm(n = NN, mean = 0, sd = Omega_Time_to_Effect))  + Delay_between_baseline_and_admin

#Max_Time_to_Regrowth, allometrically scaled 
#We simulate 1000 values, then add to each of them the delay between baseline
#assessment and administration

set.seed(12)
Theta_Max_Time_to_Regrowth_scaled = 88.5; Omega_Max_Time_to_Regrowth = 0.16 
Max_Time_to_Regrowth_scaled = Theta_Max_Time_to_Regrowth_scaled * exp(rnorm(n = NN, mean = 0, sd = Omega_Max_Time_to_Regrowth)) + Delay_between_baseline_and_admin


#hill_TtR, no variability
set.seed(13)
Theta_hill_TtR = 1.53; Omega_hill_TtR = 0.1
hill_TtR = Theta_hill_TtR * exp(rnorm(n = NN, mean = 0, sd = Omega_hill_TtR))

#Cav50_TtR
set.seed(14)
Theta_Cav50_TtR = 0.119; Omega_Cav50_TtR = 0.4
Cav50_TtR = Theta_Cav50_TtR * exp(rnorm(n = NN, mean = 0, sd = Omega_Cav50_TtR))

#Max_KL_resistant, its theta is always a 101% of the KL_sensible (proportionality)
set.seed(15)
Theta_Max_KL_resistant_medium = Theta_KL_sensible_medium * 1.01; Omega_Max_KL_resistant = 0.1
Max_KL_resistant_medium = Theta_Max_KL_resistant_medium * exp(rnorm(n = NN, mean = 0, sd = Omega_Max_KL_resistant))

#Cav50_KLres
set.seed(16)
Theta_Cav50_KLres = 0.0964; Omega_Cav50_KLres = 0.96 
Cav50_KLres = Theta_Cav50_KLres * exp(rnorm(n = NN, mean = 0, sd = Omega_Cav50_KLres))


######################
#  Combo parameters  #
######################

#Max_KL_resistant_combo, its theta is always a 32.6% of the KL_sensible (proportionality)
set.seed(17)
Theta_KL_resistant_medium_combo = Theta_KL_sensible_medium * 0.326; Omega_KL_resistant = 0.1
KL_resistant_medium_combo = Theta_KL_resistant_medium_combo * exp(rnorm(n = NN, mean = 0, sd = Omega_KL_resistant))


#Max_Time_to_Regrowth_combo, increase of 25.2% of Max_time to regrowth in monotherapy 
#We simulate 1000 values, then add to each of them the delay between baseline
#assessment and administration

set.seed(18)
Theta_Time_to_Regrowth_combo = (Theta_Max_Time_to_Regrowth_scaled * 1.252); Omega_Time_to_Regrowth_combo = 0.33
Time_to_Regrowth_combo = Theta_Time_to_Regrowth_combo * exp(rnorm(n = NN, mean = 0, sd = Omega_Time_to_Regrowth_combo)) + Delay_between_baseline_and_admin


# All the parameters in a dataframe, bearing in mind that we will
# have to change the columns every time we simulate a different condition (mono vs combo)

Parameters <- data.frame(CL, Q, V1, V2, sTYRP1, kon, koff, #PK
                         IC50, KL_sensible = KL_sensible_medium, Kill_max = Kill_max_prop_medium, 
                         Time_to_Effect = Time_to_Effect_scaled, 
                         Max_Time_to_Regrowth = Max_Time_to_Regrowth_scaled, 
                         hill_TtR, Cav50_TtR, Max_KL_resistant = Max_KL_resistant_medium, 
                         Cav50_KLres) #PD



Params_names <- names(Parameters)


#No correlation between EBEs was observed in model building, and as so,
# no correlation is included between translated params. This line checks
# which is the greatest corr in the matrix, excluding the diagonals (1),
# and the NAs that arise because of no variability in Kon/Koff. A warning
#should be thrown, because Kon/Koff have sd of 0, and max val should be "low"

#' Corr_matrix <- cor(Parameters)
#' cat(max(abs(Corr_matrix[which(Corr_matrix != 1)]), na.rm = T))

