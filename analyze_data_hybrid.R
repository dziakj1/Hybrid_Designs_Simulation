analyze_data_hybrid <- function(data_wide,
                                n_obs,
                                K,
                                print_results=TRUE) {
  # Weight and replicate;
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
  # Convert wide to long;
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
  # Fit analysis model; 
  gee_formula_proximal <- Y ~ Z1 + Z2stage2 + Z1*Z2stage2 + 
    A + A*Z1 + A*Z2stage2 + A*Z1*Z2stage2 ;
  gee_proximal <- geese(formula = gee_formula_proximal,
                        id=id, 
                        weights = known_weight, 
                        waves=wave, 
                        data=data_long_w_r,
                        corstr = "independence"); 
  #### Analysis for distal outcome ################
  # Distal outcome is assumed to be sum of proximals
  all_Y_cols <- which(substr(colnames(data_wide_w_r),1,1)=="Y");
  all_A_cols <- which(substr(colnames(data_wide_w_r),1,1)=="A");
  stage2_A_cols <- (which(colnames(data_wide_w_r)=="A1")):
    (which(colnames(data_wide_w_r)==paste("A",K+1,sep="")));
  data_wide_w_r$sum_Y <- apply(data_wide_w_r[,all_Y_cols],1,sum);
  data_wide_w_r$mean_A <- apply(data_wide_w_r[,all_A_cols],1,mean);
  data_wide_w_r$mean_A_stage2 <- apply(data_wide_w_r[,stage2_A_cols],1,mean);
  data_wide_w_r$mean_A_stage2_Z2 <- data_wide_w_r$mean_A_stage2 * 
    data_wide_w_r$Z2;
  data_wide_w_r$mean_A_stage2_Z1_Z2<- data_wide_w_r$mean_A_stage2 *  
    data_wide_w_r$Z1 * 
    data_wide_w_r$Z2;
  gee_formula_distal <- sum_Y ~ Z1 + Z2 + Z1*Z2 + 
    mean_A + mean_A*Z1 + mean_A_stage2_Z2 + mean_A_stage2_Z1_Z2 ;
  gee_distal <- geese(formula = gee_formula_distal,
                      id=id, 
                      weights = known_weight, 
                      waves=replicant, 
                      data=data_wide_w_r,
                      corstr = "independence"); 
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