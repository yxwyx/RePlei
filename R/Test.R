set.seed(123)  # for reproducibility
load("~/thesis/RePlei_Analysis/simulations/RePlei/new_sims_M5000/Cg_S1/Ce_S1/Sim_M5000_K2_Q3_p0.1_lambda0.1_alpha2_tau25e-04_CgS1_CeS1_rep1.RData")


betajk=sim_data$betajk
sjk=sim_data$sjk
C_g=sim_data$C_g
C_e=sim_data$C_e
h2=sim_data$h2

na_mask <- array(runif(5000 * 2 * 3) < 0.05, dim = c(5000, 2, 3))
betajk[na_mask] <- NA
sjk[na_mask] <- NA