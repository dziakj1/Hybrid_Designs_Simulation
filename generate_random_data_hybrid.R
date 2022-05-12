generate_random_data_hybrid <- function( n_sim,  # number of datasets to simulate;
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
                                         rho = .5) # autocorrelation between consecutive error terms in proximal outcomes;
 # Generate a simulated dataset representing a possible hybrid SMART-MRT study 
 # with hypothesized true parameters;
 { 
  library(dplyr);
  library(tidyr);
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
  return(data_wide);
}