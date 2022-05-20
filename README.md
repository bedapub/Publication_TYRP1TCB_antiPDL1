# Publication_TYRP1TCB_antiPDL1
This repository contains the R scripts to generate 1000 virtual patients, simulate their expected tumor growth profile, and evaluate it using RECIST v1.1, outputting ORR, DoR and PFS

Requires: `dplyr`, `RxODE`, `pracma` 
Whole workflow can be executed by running the script `Evaluate_Simulations.R`

`Simulate_Parameters.R` creates the virtual patients. Parameter values (both Thetas and Omegas) are hard-coded; details on how they are translated from the estimated mouse parameters are available in Supplemmentary Materials & Methods

`Run_Simulations.R` sources the `Simulate_Parameters.R` file, and uses the parameters to simulate the tumor volume changes for 400 days under 150 mg q3w of TYRP1-TCB, both alone and in combination with anti-PD-L1.

`Evaluate_Simulations.R` uses RECIST v1.1 criteria to define progression and response for the simulated tumor volume profiles, both with the full granularity of the simulation and under different tumor scanning schedules. 
It then calculates the median, 5th and 95th percentiles of the DoR and PFS, which it then prints in the console, together with the ORR, both for monotherapy, and for combination
