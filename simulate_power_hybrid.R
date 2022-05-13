simulate_power_hybrid <- function(n_sim,  # number of datasets to simulate;
                                  n_sub,  # number of subjects per simulated dataset;
                                  n_obs = 112, # number of observation occasions;
                                  K = 28, # re-randomization occasion;
                                  p_responder = .5,  # proportion of responders;
                                   # in this version we assume that it is the
                                   # same regardless of treatment arms;
                                   # The parameter values below are set by default to those used in the paper:;
                                  b0 = .25, # coefficient representing an intercept constant;
                                  b1 = -.03,  # coefficient representing half the main effect of stage-1 treatment Z1;
                                  b2 = -.03,  # coefficient representing half the effect of stage-2 treatment Z2 on nonresponders;
                                  b3 = -.03,  # coefficient representing the interaction between stage-1 and stage-2 treatment;
                                  gamma0 = -.02, # coefficient representing half the proximal main effect of intensively randomized treatment A;
                                  gamma1 = -.02, # coefficient representing interaction between A and Z1;
                                  gamma2 = -.02, # coefficient representing interaction between A and Z2;
                                  gamma3 = -.02, # coefficient representing interaction between A, Z1 and Z2;
                                  delta = -.08, # non-causal effect on proximal outcome associated with being a responder;
                                  sigma_sqd = .2, # error variance of proximal oucome;
                                  rho = .5, # autocorrelation between consecutive error terms in proximal outcomes;
                                  make_all_effects_null = FALSE, # included for convenience to override defaults and set all b's and gamma's to zero;
                                  print_results = TRUE) {
  library(dplyr);
  library(tidyr);
  start_time <- Sys.time();
  if (make_all_effects_null) {
    # A shortcut for doing Type One error simulations, 
    # to override the usual values of the parameters
    # for the randomized effects.
    b0 <- b1 <- b2 <- b3 <- 0;
    gamma0 <- gamma1 <- gamma2 <- gamma3 <- 0;
  }
  # Create empty data structures to which answers will be added;
  h0_rejected_proximal <- NULL;
  h0_rejected_distal <- NULL;
  # Start Loop;
  for (sim in 1:n_sim) {
    sim_data_wide <- generate_random_data_hybrid(n_sub = n_sub,
                                                 n_obs = n_obs,
                                                 K=K,
                                                 p_responder = p_responder,
                                                 b0 = b0,
                                                 b1 = b1,
                                                 b2 = b2,
                                                 b3 = b3,
                                                 gamma0 = gamma0,
                                                 gamma1 = gamma1,
                                                 gamma2 = gamma2, 
                                                 gamma3 = gamma3, 
                                                 delta = delta,
                                                 sigma_sqd = sigma_sqd,
                                                 rho = rho );
              # generate simulated data from specified parameters;
    results <- analyze_data_hybrid(data_wide=sim_data_wide,
                                   n_obs = n_obs,
                                   K=K,
                                   print_results = FALSE);
              # analyze the simulated data under the assumed model;
              # and return the parameter estimates and standard errors;
    # Record whether the null hypothesis was rejected or not, for
    # each coefficient in the vector of coefficients;
    h0_rejected_proximal_this_time <- summary(results$proximal)$mean$p<.05;
    names(h0_rejected_proximal_this_time) <- names(results$proximal$beta);
    h0_rejected_proximal <- bind_rows(h0_rejected_proximal,
                                      h0_rejected_proximal_this_time);
           # add the record of whether the null hypotheses were rejected
           # in this simulated dataset, to the object containing the 
           # records from any other simulations done so far.;
    h0_rejected_distal_this_time <- summary(results$distal)$mean$p<.05;
    names(h0_rejected_distal_this_time) <- names(results$distal$beta);
    h0_rejected_distal <- bind_rows(h0_rejected_distal,
                                    h0_rejected_distal_this_time);
  } # end loop over simulations;
  finish_time <- Sys.time();
  time_taken <- difftime(finish_time,start_time);
  # The code below is to print the simulated power in table form.
  # The rows are rearranged from the R default for listing model
  # terms, in order to put the SMART terms first and the MRT terms
  # afterwards.;
  row_order_proximal <- c("Z1",
                          "Z2stage2",
                          "Z1:Z2stage2",
                          "A",
                          "Z1:A",
                          "Z2stage2:A",
                          "Z1:Z2stage2:A");
  row_order_distal <-   c("Z1",
                          "Z2",
                          "Z1:Z2",
                          "mean_A",
                          "Z1:mean_A",
                          "mean_A_stage2_Z2",
                          "mean_A_stage2_Z1_Z2");
  power_proximal <- apply(h0_rejected_proximal,2,mean)[row_order_proximal];
  power_distal <- apply(h0_rejected_distal,2,mean)[row_order_distal];
  if (print_results) {
    print(cbind(proximal=power_proximal,
                distal=power_distal))
  }
  return(list(power_proximal=power_proximal, 
              power_distal=power_distal,
              time_taken=difftime(finish_time,start_time)));
}
