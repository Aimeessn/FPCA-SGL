# Simulate datasets with Model 1 (LME with a linear time trend) under 5 scenarios (as listed in Table 2) 
# n=2000 and rep=200

library(dplyr)
library(tidyverse)
library(lubridate)
library(truncnorm)
library(Hmisc)
library(faux)
library(matrixcalc)
library(MASS)
library(copula)
library(ReIns)

# Sample size = 2000
sample_size = 2000

# Repetitions = 200
rep = 200



################################################################################################################################
# Simulate coefficients for biomarkers' true trajectories Z_ki(t) for k=1,2,3,4
# 4 biomarkers in 2 groups: 1st group with relatively lower cor; 2nd group with relatively higher cor

num_bm = 4 #number of biomarkers
fixed = data.frame(intercepts = c(0.5, -0.8, 0.2, -0.5),
                   slopes = c(3, -1, 1, 1))
Sigma_low = 0.3^2*matrix(c(1, 0.2, 0.1, 0.1,
                           0.2, 1, 0.1, 0.1,
                           0.1, 0.1, 1, 0.2,
                           0.1, 0.1, 0.2, 1), nrow = 4, ncol = 4)
Sigma_high = 0.3^2*matrix(c(1, 0.8, 0.1, 0.1,
                            0.8, 1, 0.1, 0.1,
                            0.1, 0.1, 1, 0.8,
                            0.1, 0.1, 0.8, 1), nrow = 4, ncol = 4)
sim_random<- list() #a list of `rep` dfs, each df has `sample_size` rows and 2*`num_bm` columns (intercepts and slopes for each bm)
set.seed(9172)
for (r in 1:rep){
  random <- as.data.frame(cbind(mvrnorm(n = sample_size, mu = c(0,0,0,0), Sigma = Sigma_low),
                                mvrnorm(n = sample_size, mu = c(0,0,0,0), Sigma = Sigma_high)))
  colnames(random) <- c("bm1_intercepts", "bm2_intercepts", "bm1_slopes", "bm2_slopes",
                        "bm3_intercepts", "bm4_intercepts", "bm3_slopes", "bm4_slopes")
  sim_random[[r]] <- random
}
random_col <- colnames(sim_random[[1]])



################################################################################################################################
# Simulate survival times under different scenarios
# Model: hi(t) = h0(t) * exp(alpha1*Z1i(t) + alpha2*Z2i(t) + alpha3*Z3i(t) + alpha4*Z4i(t))

num_sce = 5

# Set up parameters for baseline hazard: h0(t) = gamma*t^{gamma-1}*exp(eta)
h0_gamma = 2
h0_eta = -2

# Set up functions for simulation survival outcomes using inverse transform sampling
weib_haztv <- function(t,alphas,fixed,random_i){
  h0_gamma * t^(h0_gamma-1) * exp(h0_eta) *
    exp(alphas[1]*(fixed$intercepts[1] + fixed$slopes[1]*t + random_i$bm1_intercepts + random_i$bm1_slopes*t) +
          alphas[2]*(fixed$intercepts[2] + fixed$slopes[2]*t + random_i$bm2_intercepts + random_i$bm2_slopes*t) +
          alphas[3]*(fixed$intercepts[3] + fixed$slopes[3]*t + random_i$bm3_intercepts + random_i$bm3_slopes*t) +
          alphas[4]*(fixed$intercepts[4] + fixed$slopes[4]*t + random_i$bm4_intercepts + random_i$bm4_slopes*t))
}
weib_Haztv_numeric <- function(t,alphas,fixed,random_i){
  integrate(weib_haztv, lower=0, upper=t, alphas=alphas, fixed=fixed, random_i=random_i)$value
}
weib_survtv_numeric <- function(t,alphas,fixed,random_i){
  exp(-weib_Haztv_numeric(t,alphas,fixed,random_i))
}
weib_invsurvtv_numeric <- function(p,alphas,fixed,random_i){
  uniroot(f = function(x){weib_survtv_numeric(x,alphas,fixed,random_i) - p},
          interval = c(0,80))$root
}

# Scenario 1: simulate outcome ~ bm group 1 (with low correlation)
sim_outcome_sce1 <- list()
set.seed(363)
for (r in 1:rep){
  surv_time <- numeric(sample_size)
  for(i in 1:sample_size){
    sim_unif <- runif(n=1, min=0, max=1)
    cat("r=", r, "i=", i, ": ", sim_unif, "\n")
    surv_time[i] <- try(weib_invsurvtv_numeric(p = sim_unif, alphas = c(1, 1, 0, 0),
                                               fixed = fixed, random_i = sim_random[[r]][i,]))
  }
  sim_outcome_sce1[[r]] <- data.frame(ID = seq(1, sample_size), t_to_death = surv_time)
}

# Scenario 2: simulate outcome ~ bm group 2 (with high correlation)
sim_outcome_sce2 <- list()
set.seed(4611)
for (r in 1:rep){
  surv_time <- numeric(sample_size)
  for(i in 1:sample_size){
    sim_unif <- runif(n=1, min=0, max=1)
    cat("r=", r, "i=", i, ": ", sim_unif, "\n")
    surv_time[i] <- try(weib_invsurvtv_numeric(p = sim_unif, alphas = c(0, 0, 1, 1),
                                               fixed = fixed, random_i = sim_random[[r]][i,]))
  }
  sim_outcome_sce2[[r]] <- data.frame(ID = seq(1, sample_size), t_to_death = surv_time)
}

# Scenario 3: simulate outcome ~ bm group 1 + bm group 2
sim_outcome_sce3 <- list()
set.seed(490)
for (r in 1:rep){
  surv_time <- numeric(sample_size)
  for(i in 1:sample_size){
    sim_unif <- runif(n=1, min=0, max=1)
    cat("r=", r, "i=", i, ": ", sim_unif, "\n")
    surv_time[i] <- try(weib_invsurvtv_numeric(p = sim_unif, alphas = c(1, 1, 1, 1),
                                               fixed = fixed, random_i = sim_random[[r]][i,]))
  }
  sim_outcome_sce3[[r]] <- data.frame(ID = seq(1, sample_size), t_to_death = surv_time)
}

# Scenario 4: simulate outcome randomly
sim_outcome_sce4 <- list()
set.seed(5906)
for (r in 1:rep){
  surv_time <- rweibull(n = sample_size, shape = h0_gamma, scale = 1)
  sim_outcome_sce4[[r]] <- data.frame(ID = seq(1, sample_size), t_to_death = surv_time)
}

# Resulting simulated survival outcomes: a list of `num_sce` lists (one for each simulation scenario),
# each list is a list of `rep` dfs, each df has `sample_size` rows and 3 columns (t_to_death, time, status)
for (r in 1:rep){
  if (!is.numeric(sim_outcome_sce1[[r]]$t_to_death)) {sim_outcome_sce1[[r]]$t_to_death <- as.numeric(levels(sim_outcome_sce1[[r]]$t_to_death))[sim_outcome_sce1[[r]]$t_to_death]}
  if (!is.numeric(sim_outcome_sce2[[r]]$t_to_death)) {sim_outcome_sce2[[r]]$t_to_death <- as.numeric(levels(sim_outcome_sce2[[r]]$t_to_death))[sim_outcome_sce2[[r]]$t_to_death]}
  if (!is.numeric(sim_outcome_sce3[[r]]$t_to_death)) {sim_outcome_sce3[[r]]$t_to_death <- as.numeric(levels(sim_outcome_sce3[[r]]$t_to_death))[sim_outcome_sce3[[r]]$t_to_death]}
  if (!is.numeric(sim_outcome_sce3_low[[r]]$t_to_death)) {sim_outcome_sce3_low[[r]]$t_to_death <- as.numeric(levels(sim_outcome_sce3_low[[r]]$t_to_death))[sim_outcome_sce3_low[[r]]$t_to_death]}
  if (!is.numeric(sim_outcome_sce3_high[[r]]$t_to_death)) {sim_outcome_sce3_high[[r]]$t_to_death <- as.numeric(levels(sim_outcome_sce3_high[[r]]$t_to_death))[sim_outcome_sce3_high[[r]]$t_to_death]}
}

sim_outcome <- list(sim_outcome_sce1, sim_outcome_sce2, sim_outcome_sce3,
                    sim_outcome_sce3_low, sim_outcome_sce3_high, sim_outcome_sce4)
rm(sim_outcome_sce1)
rm(sim_outcome_sce2)
rm(sim_outcome_sce3)
rm(sim_outcome_sce4)

# According to simulated survival times, choose time of censoring (EoFU) so that we observe ~10% death within FU time
#hist(sim_outcome[[1]][[1]]$t_to_death) #dist of simulated survival times, for scenario1 rep1
perc_death_obs = 0.1
t_cen_candidate <- numeric(num_sce)
for (s in 1:(num_sce-1)){
  t_cen_candidate[s] = mean(unlist(lapply(sim_outcome[[s]],
                                          function(df){quantile(df$t_to_death, probs = perc_death_obs,
                                                                na.rm = TRUE)})))
}
t_cen_candidate[1:(num_sce-1)] #0.6232951 0.6111046 0.5003066 0.3258017
t_cen = rep(0.6, num_sce)
for (s in 1:(num_sce-1)){
  for (r in 1:rep){
    sim_outcome[[s]][[r]] <- sim_outcome[[s]][[r]] %>%
      mutate(t_to_death = ifelse(is.na(t_to_death), 1000, t_to_death)) %>%
      mutate(time = pmin(t_to_death, t_cen[s]), status = as.numeric(t_to_death < t_cen[s]))
  }
}



################################################################################################################################
# Simulate true and observed biomarker trajectory data

num_measure = 10 #number of time points
time_exp_rate = 3 #s.t. simulated times (other than t=0) have mean tcen/time_exp_rate
noise_sd = 0.1 #sd for noise added to true traj to get observed traj
num_sce

sim_bm_true <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has `sample_size*num_measure` rows and 8 cols (true trajectories)
sim_bm_obs <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has `sample_size*num_measure` rows and 8 cols (observed trajectories)
set.seed(173)
for (s in 1:num_sce){
  sim_bm_true_s <- list() #a list of `rep` dfs, each df has `sample_size*num_measure` rows and 8 cols (true trajectories), one of this list for each simulation scenario
  sim_bm_obs_s <- list() #a list of `rep` dfs, each df has `sample_size*num_measure` rows and 8 cols (observed trajectories), one of this list for each simulation scenario
  for (r in 1:rep){
    # Simulate sparse and irregular time points for patients' biomarker measurements
    # with a truncated (within FU, thus truncated at tcen) Exp distribution with mean at tcen/3 (based on real data)
    sim_bm <- data.frame(ID = rep(seq(1, sample_size), each = num_measure-1),
                         Time = ReIns::rtexp(n = sample_size*(num_measure-1), rate = time_exp_rate/t_cen[s],
                                             endpoint = t_cen[s])) %>%
      arrange(ID, Time) %>%
      mutate(Obs_perID = rep(seq(2, num_measure), sample_size))%>%
      rbind(data.frame(ID = seq(1, sample_size),
                       Time = rep(0, sample_size),
                       Obs_perID = rep(1, sample_size))) %>%
      arrange(ID, Obs_perID) %>%
      mutate(Obs = seq(1, sample_size*num_measure)) %>%
      mutate(Time = round(Time, 2)) %>%
      group_by(ID) %>% distinct(Time, .keep_all=T) %>% ungroup()
    
    # Simulate biomarkers' true trajectory measurements: Z_ki(t) for k=1,2,3,4
    random_r <- sim_random[[r]][,c(1:8)]
    colnames(random_r) <- c("bm1_intercepts", "bm2_intercepts", "bm1_slopes", "bm2_slopes",
                            "bm3_intercepts", "bm4_intercepts", "bm3_slopes", "bm4_slopes") 
    bm_true <- sim_bm %>% left_join(random_r %>% mutate(ID = seq(1, sample_size)), by = "ID") %>%
      mutate(bm1 = fixed$intercepts[1] + fixed$slopes[1]*Time + bm1_intercepts + bm1_slopes*Time,
             bm2 = fixed$intercepts[2] + fixed$slopes[2]*Time + bm2_intercepts + bm2_slopes*Time,
             bm3 = fixed$intercepts[3] + fixed$slopes[3]*Time + bm3_intercepts + bm3_slopes*Time,
             bm4 = fixed$intercepts[4] + fixed$slopes[4]*Time + bm4_intercepts + bm4_slopes*Time)
    bm_true_sd = mean(c(sd(bm_true$bm1), sd(bm_true$bm2), sd(bm_true$bm3), sd(bm_true$bm4)))
    sim_bm_true_s[[r]] <- sim_bm %>% cbind(bm_true %>% dplyr::select(bm1, bm2, bm3, bm4))
    
    # Simulate biomarkers' observed trajectories measurements (true + noise + missingness)
    noise = matrix(rnorm(n = nrow(sim_bm)*num_bm, mean = 0, sd = noise_sd),
                   nrow = nrow(sim_bm), ncol = num_bm)
    sim_bm_obs_s_temp <- sim_bm %>% cbind(bm_true %>% dplyr::select(bm1, bm2, bm3, bm4) + noise)
    
    # Simulate missingness due to death and censoring from discharge
    gamma_shape = 2
    gamma_scale = t_cen[s]/8*3
    if (s<=4){
      sim_bm_obs_s_temp <- sim_bm_obs_s_temp %>%
        left_join(sim_outcome[[s]][[r]] %>%
                    mutate(t_hosp = rgamma(n = sample_size, shape = gamma_shape, scale = gamma_scale)), by = "ID") %>%
        mutate(bm1_missing = case_when(status==1 & Time<=t_to_death ~ bm1,
                                       status==1 & Time>t_to_death ~ NA_real_,
                                       status==0 & Time<=t_hosp ~ bm1,
                                       T ~ NA_real_),
               bm2_missing = case_when(status==1 & Time<=t_to_death ~ bm2,
                                       status==1 & Time>t_to_death ~ NA_real_,
                                       status==0 & Time<=t_hosp ~ bm2,
                                       T ~ NA_real_),
               bm3_missing = case_when(status==1 & Time<=t_to_death ~ bm3,
                                       status==1 & Time>t_to_death ~ NA_real_,
                                       status==0 & Time<=t_hosp ~ bm3,
                                       T ~ NA_real_),
               bm4_missing = case_when(status==1 & Time<=t_to_death ~ bm4,
                                       status==1 & Time>t_to_death ~ NA_real_,
                                       status==0 & Time<=t_hosp ~ bm4,
                                       T ~ NA_real_))
    }
    else if (s==5){
      sim_bm_obs_s_temp <- sim_bm_obs_s_temp %>%
        mutate(bm1_missing = bm1,
               bm2_missing = bm2,
               bm3_missing = bm3,
               bm4_missing = bm4)
    }
    
    sim_bm_obs_s[[r]] <- sim_bm_obs_s_temp %>%
      dplyr::select(ID, Time, Obs_perID, Obs, bm1_missing, bm2_missing, bm3_missing, bm4_missing) %>%
      rename(bm1 = bm1_missing, bm2 = bm2_missing, bm3 = bm3_missing, bm4 = bm4_missing)
  }
  sim_bm_true <- append(sim_bm_true, list(sim_bm_true_s))
  sim_bm_obs <- append(sim_bm_obs, list(sim_bm_obs_s))
}

sapply(sim_bm_obs[[1]][[1]], function(x) sum(is.na(x)/nrow(sim_bm_obs[[1]][[1]])))
sapply(sim_bm_obs[[2]][[1]], function(x) sum(is.na(x)/nrow(sim_bm_obs[[2]][[1]])))
sapply(sim_bm_obs[[3]][[1]], function(x) sum(is.na(x)/nrow(sim_bm_obs[[3]][[1]])))
sapply(sim_bm_obs[[4]][[1]], function(x) sum(is.na(x)/nrow(sim_bm_obs[[4]][[1]])))
sapply(sim_bm_obs[[5]][[1]], function(x) sum(is.na(x)/nrow(sim_bm_obs[[5]][[1]])))



save(sim_bm_obs, sim_outcome, file = "Simulation_survival_datasets.RData")

