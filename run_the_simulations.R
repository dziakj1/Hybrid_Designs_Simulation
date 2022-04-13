rm(list = ls());
library(dplyr);
library(tidyr);
library(geepack);
set.seed(16802);
source("generate_random_data_hybrid.R");
source("analyze_data_hybrid.R");
source("simulate_power_hybrid.R");
n_sim <- 2000;
scenario <- 0;
overall_start_time <- Sys.time();
proximal_T1E <- NULL;
distal_T1E <- NULL;
proximal_power <- NULL; 
distal_power <- NULL;
for (n_sub in c(100,150,200)) {
  for (p_responder in c(.6,.5,.4)) {
    scenario <- scenario + 1;
    print("Simulating Type One error when all effects are null");
    simulation_results <- simulate_power_hybrid(n_sim=n_sim,
                                                p_responder=p_responder,
                                                n_sub=n_sub,
                                                make_all_effects_null = TRUE); 
    proximal_T1E <- rbind(proximal_T1E,
                          c(n_sub=n_sub,
                            p_responder=p_responder,
                            simulation_results$power_proximal));
    distal_T1E <- rbind(distal_T1E,
                        c(n_sub=n_sub,
                          p_responder=p_responder,
                          simulation_results$power_distal));
    print("Simulating power when all effects are non-null");
    simulation_results <- simulate_power_hybrid(n_sim=n_sim,
                                                p_responder=p_responder,
                                                n_sub=n_sub,
                                                make_all_effects_null = FALSE); 
    proximal_power <- rbind(proximal_power,
                            c(n_sub=n_sub,
                              p_responder=p_responder,
                              simulation_results$power_proximal));
    distal_power <- rbind(distal_power,
                          c(n_sub=n_sub,
                            p_responder=p_responder,
                            simulation_results$power_distal));
    save.image("okay_so_far.rdata");
  }
}
overall_finish_time <- Sys.time();
save.image("did-simulations.rdata");
print(round(proximal_T1E,4));
print(round(distal_T1E,4));
print(round(proximal_power,4));
print(round(distal_power,4));
