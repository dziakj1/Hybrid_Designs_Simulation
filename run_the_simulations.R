set.seed(16802); # set arbitrary seed for generating random data;
source("generate_random_data_hybrid.R");
source("analyze_data_hybrid.R");
source("simulate_power_hybrid.R");
n_sim <- 2500; # do 2500 simulations;
overall_start_time <- Sys.time(); # record starting time 
                                  #  just to see how long
                                  #  the scenario takes.;
# Initialize data structures that will later hold the results:
proximal_T1E <- NULL;
distal_T1E <- NULL;
proximal_power <- NULL; 
distal_power <- NULL;
# Begin loop over simulation scenarios.  We will consider
# three different sample sizes and three different 
# response probabilities (proportions of responders, 
# set for simplicity to be the same for each treatment arm).;
for (n_sub in c(100,150,200)) {
  for (p_responder in c(.6,.5,.4)) {
    # Do power simulation first with all treatment effects
    # set to zero -- this gives estimates for Type One error
    # under such a scenario.;
    print("Simulating Type One error when all effects are null");
    simulation_results <- simulate_power_hybrid(n_sim=n_sim,
                                                p_responder=p_responder,
                                                n_sub=n_sub,
                                                make_all_effects_null = TRUE); 
    # Record the results of the Type I error simulation by
    # appending the proximal and distal type I error estimates
    # for each variable as a new row in the appropriate data frame
    # of results.  Also include the number of subjects and proportion
    # of responders at the beginning of the row, so that it can be
    # read as a results table. There is a separate 'table' for 
    # proximal and for distal Type I error.;
    proximal_T1E <- rbind(proximal_T1E,
                          c(n_sub=n_sub,
                            p_responder=p_responder,
                            simulation_results$power_proximal));
    distal_T1E <- rbind(distal_T1E,
                        c(n_sub=n_sub,
                          p_responder=p_responder,
                          simulation_results$power_distal));
        # records the distal Type I error rate; 
    # Now simulate power with the nonzero default parameter values.;
    print("Simulating power when all effects are non-null");
    simulation_results <- simulate_power_hybrid(n_sim=n_sim,
                                                p_responder=p_responder,
                                                n_sub=n_sub,
                                                make_all_effects_null = FALSE); 
    proximal_power <- rbind(proximal_power,
                            c(n_sub=n_sub,
                              p_responder=p_responder,
                              simulation_results$power_proximal));
        # records the proximal power results;
    distal_power <- rbind(distal_power,
                          c(n_sub=n_sub,
                            p_responder=p_responder,
                            simulation_results$power_distal));
        # records the distal power results;
    save.image("okay_so_far.rdata");
  }
}
# Save and print answers:;
overall_finish_time <- Sys.time();
save.image("did-simulations.rdata");
print(round(proximal_T1E,4));
print(round(distal_T1E,4));
print(round(proximal_power,4));
print(round(distal_power,4));
