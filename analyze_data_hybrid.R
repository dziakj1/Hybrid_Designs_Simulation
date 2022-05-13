analyze_data_hybrid <- function(data_wide, # the dataset in wide format;
                                n_obs,  # number of observations total (e.g., number of days if observations are daily);
                                K,  # time (e.g., day) of second randomization;
                                print_results=TRUE) {
  # Do the analysis on a simulated dataset as part of the power simulation.;
  # Might also be useful for real data in the same format. 
  library(dplyr);
  library(tidyr);
  library(geepack);
  # Weight and replicate;
  # Responders don't actually get randomized to the stage-2 factor Z2,
  # so their data is replicated into two copies, one with each possible 
  # value of Z2 which they could have received had they been a responder. 
  # To counterbalance this, nonresponders get double weight (4 instead
  # of 2) in the final analysis. This can be interpreted as inverse 
  # propensity weighting, although a special kind used in a randomized
  # trial. See, e.g.,
  # Daniel Almirall, Inbal Nahum-Shani, Nancy E. Sherwood, Susan A. Murphy, 
  # Introduction to SMART designs for the development of adaptive 
  # interventions: with application to weight loss research, 
  # Translational Behavioral Medicine, Volume 4, Issue 3, September 2014, 
  # pages 260-274, https://doi.org/10.1007/s13142-014-0265-0
  # or
  # Xi Lu, Inbal Nahum-Shani, Connie Kasari, Kevin G. Lynch, David W. Oslin, 
  # William E. Pelham, Gregory Fabiano, Daniel Almirall. Comparing dynamic 
  # treatment regimes using repeated-measures outcomes: modeling 
  # considerations in SMART studies. Statistics in Medicine, Volume 35, 
  # Issue10, May 2016, Pages 1595-1615, 
  # https://onlinelibrary.wiley.com/doi/full/10.1002/sim.6819
  rows_to_replicate <- data_wide %>% filter(R==1);
  rows_to_replicate$known_weight <- 2;
  positive_pseudodata <- rows_to_replicate;
  positive_pseudodata$Z2 <- +1;
  positive_pseudodata$replicant <- 1L;
  negative_pseudodata <- rows_to_replicate;
  negative_pseudodata$Z2 <- -1;
  negative_pseudodata$replicant <- 2L;
  rows_not_to_replicate <- data_wide %>% filter(R==0);
  rows_not_to_replicate$replicant <- 1;
  rows_not_to_replicate$known_weight <- 4;
  data_wide_w_r <- bind_rows(rows_not_to_replicate,
                             positive_pseudodata,
                             negative_pseudodata);
  data_wide_w_r <- data_wide_w_r %>% 
    arrange(id, replicant);
  #### Analysis for proximal outcome ################
  # Convert wide dataset to long dataset (i.e., now
  # have rows for each observation within each person,
  # instead of just one row per person) using 
  # tidyverse tricks to make it simpler;
  data_long_w_r_A <- data_wide_w_r %>% 
    select("id","replicant",starts_with("A")) %>%
    pivot_longer(cols=starts_with("A"),
                 names_prefix = "A",
                 names_to = "time", 
                 values_to = "A");
  data_long_w_r_Y <- data_wide_w_r %>% 
    select("id","replicant",starts_with("Y")) %>%
    pivot_longer(cols=starts_with("Y"),
                 names_prefix = "Y",
                 names_to = "time", 
                 values_to = "Y");
  data_long_w_r_other <- data_wide_w_r %>% select("id", 
                                                  "replicant",
                                                  "known_weight",
                                                  "Z1",
                                                  "R",
                                                  "Z2");
  data_long_w_r <- inner_join(data_long_w_r_other,
                              data_long_w_r_A, 
                              by=c("id","replicant")) %>% 
    inner_join(data_long_w_r_Y,by=c("id","replicant","time"))  %>% 
    mutate(wave = as.integer(time) + n_obs*(replicant-1)) %>%
    mutate(stage2 = 1*(time>K)) %>%
    mutate(Z2stage2=Z2*stage2);
  # Fit analysis model using generalized estimating software
  # in the geepack package. See Liang and Zeger (1986, Biometrika).
  # geepack is available from CRAN, contributed by 
  # Soren Hojsgaard, Ulrich Halekoh, Jun Yan, and Claus Ekstrom, at
  # https://cran.r-project.org/web/packages/geepack/index.html ;
  gee_formula_proximal <- Y ~ Z1 + Z2stage2 + Z1*Z2stage2 + 
    A + A*Z1 + A*Z2stage2 + A*Z1*Z2stage2 ;
  gee_proximal <- geese(formula = gee_formula_proximal,
                        id=id, 
                        weights = known_weight, 
                        waves=wave, 
                        data=data_long_w_r,
                        corstr = "independence"); 
  #### Analysis for distal outcome ################
  # Distal outcome is assumed to be sum of proximal outcomes. This
  # will not be true in all studies but is assumed to be true in the
  # current example.;
  all_Y_cols <- which(substr(colnames(data_wide_w_r),1,1)=="Y");
  all_A_cols <- which(substr(colnames(data_wide_w_r),1,1)=="A");
  stage2_A_cols <- (which(colnames(data_wide_w_r)=="A1")):
    (which(colnames(data_wide_w_r)==paste("A",K+1,sep="")));
  data_wide_w_r$sum_Y <- apply(data_wide_w_r[,all_Y_cols],1,sum);
      # The sum of Y values for each individual is computed, to use as the 
       # distal outcome.;
  data_wide_w_r$mean_A <- apply(data_wide_w_r[,all_A_cols],1,mean);
  data_wide_w_r$mean_A_stage2 <- apply(data_wide_w_r[,stage2_A_cols],1,mean);
  data_wide_w_r$mean_A_stage2_Z2 <- data_wide_w_r$mean_A_stage2 * 
    data_wide_w_r$Z2;
  data_wide_w_r$mean_A_stage2_Z1_Z2<- data_wide_w_r$mean_A_stage2 *  
    data_wide_w_r$Z1 * 
    data_wide_w_r$Z2;
      # The mean value of the multiply randomized factor A is used as a
      # covariate;  it may be easier to interpret than the sum.;
  gee_formula_distal <- sum_Y ~ Z1 + Z2 + Z1*Z2 + 
    mean_A + mean_A*Z1 + mean_A_stage2_Z2 + mean_A_stage2_Z1_Z2 ;
  # The distal is related to its appropriate covariates via a separate model:
  gee_distal <- geese(formula = gee_formula_distal,
                      id=id, 
                      weights = known_weight, 
                      waves=replicant, 
                      data=data_wide_w_r,
                      corstr = "independence"); 
  # Now return the results:
  results <- list(proximal=gee_proximal,
                  distal=gee_distal);
  if (print_results) {
    print("Proximal model results:")
    print(summary(results$proximal)$mean);
    print("Distal model results:")
    print(summary(results$distal)$mean);
  }
  return(results);
}