rm(list = ls());
library(dplyr);
library(tidyr);
library(geepack);
set.seed(16801); 
n_sim <- 1000;
simulate_power <- function(active,
                           n_sub, 
                           n_sim,p_responder=p_responder) {
  start_time <- Sys.time(); 
  K <- 28; # re-randomization occasion;
  n_obs <- 112;
  b0 <- .25;
  if (active) {
  b1 <- b2 <- b3 <- -.03;
  gamma0 <- gamma1 <- gamma2 <- gamma3 <- -.02; 
  } else {
    b1 <- b2 <- b3 <- 0;
    gamma0 <- gamma1 <- gamma2 <- gamma3 <- 0; 
  }
  delta <- -.08; 
  sigma_sqd <- .2;
  rho <- .5; 
  # Create data structure to hold answers
  h0_rejected_proximal <- NULL;
  h0_rejected_distal <- NULL;
  # Start Loop;
  for (sim in 1:n_sim) {
    # Generate data for Stage 1;
    data_wide <- data.frame(id=1:n_sub,
                            replicant=NA,
                            known_weight=NA,
                            Z1=sample(c(1,-1),
                                      n_sub,
                                      replace=TRUE));
    n_obs_stage1 <- K-1;
    A_stage1_wide <- matrix(sample(c(1,-1),
                                   n_sub*n_obs_stage1,
                                   replace=TRUE),n_sub,n_obs_stage1);
    Z1_stage1_wide <- data_wide$Z1 %o% rep(1,n_obs_stage1);
    mu_prox_stage1 <- b0 + b1*Z1_stage1_wide + gamma0*A_stage1_wide;
    resids_stage1 <- matrix(NA,
                            nrow(mu_prox_stage1),
                            ncol(mu_prox_stage1));
    resids_stage1[,1] <- rnorm(n_sub,
                               mean=0,
                               sd=sqrt(sigma_sqd));
    for (j in 2:ncol(mu_prox_stage1)) {
      resids_stage1[,j] <- rho*resids_stage1[,j-1] +
        sqrt(1-rho^2)*rnorm(n_sub,
                            mean=0,
                            sd=sqrt(sigma_sqd));
    }
    Y_prox_stage1 <- mu_prox_stage1 + resids_stage1;
    # Generate response status;
    data_wide$R <- NA;
    data_wide$R[which(data_wide$Z1==+1)] <- rbinom(sum(data_wide$Z1==+1),
                                                   1,
                                                   p_responder);
    data_wide$R[which(data_wide$Z1==-1)] <- rbinom(sum(data_wide$Z1==-1),
                                                   1,
                                                   p_responder);
    # Generate data for Stage 2;
    n_obs_stage2 <- n_obs-K+1;
    A_stage2_wide <- matrix(sample(c(1,-1),
                                   n_sub*n_obs_stage2,
                                   replace=TRUE),
                            n_sub,
                            n_obs_stage2);
    data_wide$Z2 <- 0;
    nonresponders_id <- which(data_wide$R==0);
    data_wide$Z2[nonresponders_id] <- sample(c(1,-1),
                                             length(nonresponders_id),replace=TRUE);
    Z1_stage2_wide <- data_wide$Z1 %o% rep(1,n_obs_stage2);
    Z2_stage2_wide <- data_wide$Z2 %o% rep(1,n_obs_stage2);
    mu_prox_stage2 <- b0 + b1*Z1_stage2_wide + b2*Z2_stage2_wide+ 
      b3*Z1_stage2_wide*Z2_stage2_wide +
      gamma0*A_stage2_wide + 
      gamma1*A_stage2_wide*Z1_stage2_wide + 
      gamma2*A_stage2_wide*Z2_stage2_wide +
      gamma3*A_stage2_wide*Z1_stage2_wide*Z2_stage2_wide +
      delta*(data_wide$R-(p_responder-1))  ;
    resids_stage2 <- matrix(NA,
                            nrow(mu_prox_stage2),
                            ncol(mu_prox_stage2));
    resids_stage2[,1] <- rnorm(n_sub,
                               mean=0,
                               sd=sqrt(sigma_sqd));
    for (j in 2:ncol(mu_prox_stage2)) {
      resids_stage2[,j] <- rho*resids_stage2[,j-1] +
        sqrt(1-rho^2)*rnorm(n_sub,
                            mean=0,
                            sd=sqrt(sigma_sqd));
    }
    Y_prox_stage2 <- mu_prox_stage2 + resids_stage2;
    # Assemble wide dataset;
    A_matrix <- data.frame(id=data_wide$id,
                           A=cbind(A_stage1_wide, 
                                   A_stage2_wide));
    colnames(A_matrix) <- c("id",paste("A",1:n_obs,sep=""));
    Y_matrix <- data.frame(id=data_wide$id,
                           Y=cbind(Y_prox_stage1, 
                                   Y_prox_stage2));
    colnames(Y_matrix) <- c("id",paste("Y",1:n_obs,sep=""));
    data_wide <- left_join(data_wide, 
                           A_matrix, 
                           by="id");
    data_wide <- left_join(data_wide, 
                           Y_matrix, 
                           by="id");
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
    # Record results for power analysis;
    h0_rejected_proximal_this_time <- summary(gee_proximal)$mean$p<.05;
    names(h0_rejected_proximal_this_time) <- names(gee_proximal$beta);
    h0_rejected_proximal <- bind_rows(h0_rejected_proximal,
                                      h0_rejected_proximal_this_time);
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
    h0_rejected_distal_this_time <- summary(gee_distal)$mean$p<.05;
    names(h0_rejected_distal_this_time) <- names(gee_distal$beta);
    h0_rejected_distal <- bind_rows(h0_rejected_distal,
                                    h0_rejected_distal_this_time);
  } # end loop over simulations
  finish_time <- Sys.time();
  time_taken <- difftime(finish_time,start_time);
  power_proximal <- apply(h0_rejected_proximal,2,mean);
  power_distal <- apply(h0_rejected_distal,2,mean);
  save.image(paste("ran-",
                   active,"-",
                   n_sub, "-",
                   n_sim, "-",
                   10*p_responder,
                   ".rdata",sep=""));
  return(list(active=active,
              n_sub=n_sub, 
              n_sim=n_sim, 
              p_responder=p_responder,
              time_taken=time_taken,
              power_proximal=power_proximal, 
              power_distal=power_distal));
}

scenario <- 0;
overall_start_time <- Sys.time();
answers_null <- NULL;
answers_active <- NULL;
proximal_T1E <- NULL;
distal_T1E <- NULL;
proximal_power <- NULL; 
distal_power <- NULL;
  for (n_sub in c(100,150,200)) {
    for (p_responder in c(.6,.5,.4)) {
      scenario <- scenario + 1;
      answers_null[[scenario]] <- simulate_power(active=FALSE,
                                       n_sub=n_sub, 
                                       n_sim=n_sim, 
                                       p_responder=p_responder);
      proximal_T1E <- cbind(proximal_T1E,
                            answers_null[[scenario]]$power_proximal);
      distal_T1E <- cbind(distal_T1E,
                          answers_null[[scenario]]$power_distal);
      answers_active[[scenario]] <- simulate_power(active=TRUE,
                                             n_sub=n_sub, 
                                             n_sim=n_sim, 
                                             p_responder=p_responder);
      proximal_power <- cbind(proximal_power,
                             answers_active[[scenario]]$power_proximal); 
      distal_power <- cbind(distal_power,
                             answers_active[[scenario]]$power_distal);
      save.image("oksofar-SimulationNormal4.rdata");
    }
}
overall_finish_time <- Sys.time();
cat(paste("Did",n_sim,"simulations."));

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
print(proximal_T1E[row_order_proximal,]);
print(proximal_power[row_order_proximal,]);
print(distal_T1E[row_order_distal,]);
print(distal_power[row_order_distal,]);
save.image("RanSimulationNormal4.rdata");
print(difftime(overall_finish_time,
               overall_start_time));
