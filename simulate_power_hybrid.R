simulate_power_hybrid <- function(n_sim,
                                  n_sub,
                                  n_obs = 112, # number of observation occasions;
                                  K = 28, # re-randomization occasion;
                                  p_responder = .5,
                                  b0 = .25,
                                  b1 = -.03,
                                  b2 = -.03,
                                  b3 = -.03,
                                  gamma0 = -.02,
                                  gamma1 = -.02,
                                  gamma2 = -.02, 
                                  gamma3 = -.02, 
                                  delta = -.08,
                                  sigma_sqd = .2,
                                  rho = .5,
                                  make_all_effects_null = FALSE,
                                  print_results = TRUE) {
  library(dplyr);
  library(tidyr);
  start_time <- Sys.time();
  if (make_all_effects_null) {
    # A shortcut for doing Type One error simulations, 
    # to override the usual values of the parameters
    # for the randomized effects.
    b0 <- b1 <- b2 <- b3<- 0;
    gamma0 <- gamma1 <- gamma2 <- gamma3 <- 0;
  }
  # Create data structure to hold answers
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
    results <- analyze_data_hybrid(data_wide=sim_data_wide,
                                   n_obs = n_obs,
                                   K=K,
                                   print_results = FALSE);
    # Record results for power analysis;
    h0_rejected_proximal_this_time <- summary(results$proximal)$mean$p<.05;
    names(h0_rejected_proximal_this_time) <- names(results$proximal$beta);
    h0_rejected_proximal <- bind_rows(h0_rejected_proximal,
                                      h0_rejected_proximal_this_time);
    h0_rejected_distal_this_time <- summary(results$distal)$mean$p<.05;
    names(h0_rejected_distal_this_time) <- names(results$distal$beta);
    h0_rejected_distal <- bind_rows(h0_rejected_distal,
                                    h0_rejected_distal_this_time);
  } # end loop over simulations
  finish_time <- Sys.time();
  time_taken <- difftime(finish_time,start_time);
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
