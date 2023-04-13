# Run the comparator method using baseline measurements in SGL on simulated data (3 models and 5 scenarios as listed in Table 2)
# Calculate TPR and FPR

library(dplyr)
library(tidyverse)
library(SGL)
library(parallel)
library(foreach)

load("Simulation_survival_datasets.RData")
load("Simulation_survival_quadratic_datasets.RData")
load("Simulation_survival_splines_datasets.RData")

sample_size = 2000
rep = 200 #200 repetitions
num_sce = 5

# Set up cluster with multiple cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
clusterEvalQ(my.cluster,{
  library(dplyr)
  library(SGL)
})
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()



################################################################################################################################
# Extract baseline measurement for each biomarker for each patient

bm_group = data.frame(bm = c("bm1", "bm2", "bm3", "bm4"), group = c(1,1,2,2))
baseline_func <- function(vec){
  return(first(na.omit(vec)))
}

# Data simulated with Model 1
sim_baseline <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has `sample_size` rows and 5 cols (ID+4bm)
sim_missing <- list()
for (s in 1:num_sce){
  sim_baseline_s <- list() #a list of `rep` dfs, each df has `sample_size` rows and 13 cols, one of this list for each simulation scenario
  sim_missing_s <- numeric(rep)
  for (r in 1:rep){
    bm_baseline <- sim_bm_obs[[s]][[r]] %>% group_by(ID) %>% arrange(Time, .by_group = TRUE) %>%
      dplyr::select(ID, all_of(bm_group$bm)) %>% group_by(ID) %>% summarise_all(baseline_func)
    sim_missing_s[r] <- length(table(is.na(bm_baseline)))
    sim_baseline_s[[r]] <- bm_baseline
  }
  sim_baseline <- append(sim_baseline, list(sim_baseline_s))
  sim_missing[[s]] <- sim_missing_s
}

# Data simulated with Model 2
sim_baseline_quad <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has `sample_size` rows and 5 cols (ID+4bm)
sim_missing_quad <- list()
for (s in 1:num_sce){
  sim_baseline_quad_s <- list() #a list of `rep` dfs, each df has `sample_size` rows and 13 cols, one of this list for each simulation scenario
  sim_missing_quad_s <- numeric(rep)
  for (r in 1:rep){
    bm_baseline_quad <- sim_bm_obs_quad[[s]][[r]] %>% group_by(ID) %>% arrange(Time, .by_group = TRUE) %>%
      dplyr::select(ID, all_of(bm_group$bm)) %>% group_by(ID) %>% summarise_all(baseline_func)
    sim_missing_quad_s[r] <- length(table(is.na(bm_baseline_quad)))
    sim_baseline_quad_s[[r]] <- bm_baseline_quad
  }
  sim_baseline_quad <- append(sim_baseline_quad, list(sim_baseline_quad_s))
  sim_missing_quad[[s]] <- sim_missing_quad_s
}

# Data simulated with Model 3
sim_baseline_splines <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has `sample_size` rows and 5 cols (ID+4bm)
sim_missing_splines <- list()
for (s in 1:num_sce){
  sim_baseline_splines_s <- list() #a list of `rep` dfs, each df has `sample_size` rows and 13 cols, one of this list for each simulation scenario
  sim_missing_splines_s <- numeric(rep)
  for (r in 1:rep){
    bm_baseline_splines <- sim_bm_obs_splines[[s]][[r]] %>% group_by(ID) %>% arrange(Time, .by_group = TRUE) %>%
      dplyr::select(ID, all_of(bm_group$bm)) %>% group_by(ID) %>% summarise_all(baseline_func)
    sim_missing_splines_s[r] <- length(table(is.na(bm_baseline_splines)))
    sim_baseline_splines_s[[r]] <- bm_baseline_splines
  }
  sim_baseline_splines <- append(sim_baseline_splines, list(sim_baseline_splines_s))
  sim_missing_splines[[s]] <- sim_missing_splines_s
}



################################################################################################################################
# Conduct SGL on the baseline measurements, with pre-picked alpha=0.05 (~GL), and pick lambda using CV

# Main function to perform SGL with pre-picked alpha=0.05 (~GL), and pick lambda using CV
# The range of 50 potential lambda values is defined as:
#   upper bound being the smallest lambda that shrink everything to 0, lower bound
#   lower bound is upper bound / 100
# lambda_min: the lambda that gives the minimum CV lldiff
# lambda_1se: the largest lambda (more shrinkage) that gives the CV lldiff within 1 standard error of the minimum
# Arguments: a dataframe "outcome_df" for outcome data simulated under different scenario
#            a dataframe "bm_df" for simulated biomarkers data (baseline measurements)
# Output: a list of `num_sce` lists (one for each simulation scenario),
#         each list is a list of `rep` dfs,
#         each df has p+2 (p is the number of covariates, not including intercept, the other 2 is one for
#         lambda index and one for lambda) rows and 3 columns (variable, bm_group, beta_best)
SGL_manual_func <- function(outcome_df, bm_df){
  # Combine biomarker, demographics, and outcome data
  df <- outcome_df %>% dplyr::select(ID, time, status) %>%
    left_join(bm_df, by = "ID")
  
  # Standardize biomarker data
  df[c(4:7)] <- scale(df[c(4:7)])
  
  # Prepare for SGL
  SGL_x <- model.matrix(terms( ~ bm1 + bm2 + bm3 + bm4, keep.order = TRUE),
                        df %>% dplyr::select(-ID, -time, -status))
  SGL_list <- list(x = SGL_x[,-1], time = df$time, status = df$status)
  SGL_group <- c(rep(1,2),rep(2,2)) 
  
  # Performm SGL with lambdas = NULL to get the smallest lambda that shrink everything to 0
  set.seed(483)
  SGL_temp <- cvSGL(SGL_list, index = SGL_group, type = "cox", maxit = 1000, thresh = 0.001,
                    min.frac = 0.05, nlam = 2, gamma = 0.8, nfold = 10, standardize = FALSE,
                    verbose = FALSE, step = 1, reset = 10, alpha = 0.05, foldid = NULL,
                    lambdas = NULL)
  lambda_largest <- SGL_temp$lambdas[1]
  
  # Performm SGL with a sequence of 50 lambdas, with the upper bound being the smallest lambda that shrink everything to 0s
  set.seed(483)
  SGL <- cvSGL(SGL_list, index = SGL_group, type = "cox", maxit = 1000, thresh = 0.001,
               min.frac = 0.05, nlam = 50, gamma = 0.8, nfold = 10, standardize = FALSE,
               verbose = FALSE, step = 1, reset = 10, alpha = 0.05, foldid = NULL,
               lambdas = seq(lambda_largest/100, lambda_largest*1.001, length.out = 50))
  #plot(SGL$lambdas, SGL$lldiff)
  beta_all <- SGL[["fit"]][["beta"]]
  
  # Find lambda_min and lambda_1se
  lambda_min_index <- which.min(SGL$lldiff)
  lldiff_sd = sd(SGL$lldiff)
  outliers = SGL$lldiff > min(SGL$lldiff)+2*lldiff_sd #removing outliers before calculating sd
  if (sum(outliers) <= length(SGL$lambdas)/10){ #potential outliers defined as those>min+2sd, remove if number of potential outliers is <=1/10 of total number of lambdas
    lldiff = SGL$lldiff[-which(outliers)]
    lldiff_sd = sd(lldiff)
  }
  temp = which(SGL$lldiff>=min(SGL$lldiff)+lldiff_sd)
  if (sum(temp>lambda_min_index) == 0) { #all lambda larger than lambda_min give CV lldiff within 1sd, thus lambda_1se is the largest lambda
    lambda_1se_index = length(SGL$lambdas)
  } else { #lambda_1se is the smallest lambda larger than lambda_min and gives CV lldiff within 1sd
    lambda_1se_index = min(temp[temp>lambda_min_index])
  }
  
  # If best lambda is not the largst (most shrinkage), then manually select lambda_1se
  SGL_index <- lambda_min_index
  if (SGL_index != 1){SGL_index = lambda_1se_index} 
  
  # Clean final model and save beta coeffcients and CV-picked lambda and alpha
  variable = c(colnames(df)[4:7], "index", "lambda")
  bm_group = c(rep("Group 1",2), rep("Group 2",2), "index", "lambda")
  result <- data.frame(beta_best = SGL$fit$beta[,SGL_index]) %>%
    rbind(c(SGL_index)) %>%
    rbind(c(SGL$lambdas[SGL_index])) %>%
    mutate(variable = variable, bm_group = bm_group) %>%
    dplyr::select(variable, bm_group, everything())
  
  return(result)
}

# Function to clean results from SGL_manual_func()
# Arguments: a list "result_list" straight from SGL_manual_func()
# Output: a list of 2 dfs, "beta_clean" about estimated beta coefficients with the best lambda,
#         and "index_lambda_lldiff" for the number of times the best lambda is the largest lambda and the mean lambda
clean_func <- function(result_list){
  index = 50
  beta_clean <- result_list[[1]][[1]] %>%
    filter(bm_group %in% c("Group 1", "Group 2")) %>% dplyr::select(variable, bm_group)
  index_lambda <- data.frame(Scenario1 = rep(NA, 2), Scenario2 = rep(NA, 2), Scenario3 = rep(NA, 2),
                             Scenario3low = rep(NA, 2), Scenario3high = rep(NA, 2), Scenario4 = rep(NA, 2),
                             Scenario4null = rep(NA, 2))
  
  for (s in 1:num_sce){
    beta = matrix(nrow = nrow(beta_clean), ncol = rep)
    beta_non0 = rep(0, nrow(beta_clean))
    index_numlargest = 0
    lambda_s = rep(NA, rep)
    for (r in 1:rep){
      beta[,r] <- result_list[[s]][[r]]$beta_best[-c(5:6)]
      if (result_list[[s]][[r]]$beta_best[5] == index) {index_numlargest = index_numlargest + 1}
      lambda_s[r] = result_list[[s]][[r]]$beta_best[6]
      for (i in 1:length(beta_non0)){
        if (result_list[[s]][[r]]$beta_best[i] != 0) {beta_non0[i] = beta_non0[i] + 1}
      }
    }
    beta_clean <- beta_clean %>%
      cbind(paste(format(round(apply(beta, 1, mean),4), nsmall = 4), " (",
                  format(round(apply(beta, 1, sd),4), nsmall = 4), ")", sep="")) %>%
      cbind(paste(format(round(beta_non0/rep*100, 1), nsmall=1), "%", sep=""))
    index_lambda[1,s] <- index_numlargest
    index_lambda[2,s] <- round(mean(lambda_s),4)
  }
  colnames(beta_clean)[3:12] <- c("Scenario1_beta", "Scenario1_power_size",
                                  "Scenario2_beta", "Scenario2_power_size",
                                  "Scenario3_beta", "Scenario3_power_size",
                                  "Scenario4_beta", "Scenario4_power_size",
                                  "Scenario5_beta", "Scenario5_power_size")
  rownames(index_lambda) <- c("index_numlargest", "lambda")
  
  return(list(beta_clean = beta_clean, index_lambda = index_lambda))
}

# Run SGL for each scenario and each rep on data simulated with Model 1
SGL_result <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has 5 rows and 3 cols for the CV-picked SGL beta coefficients
for (s in 1:num_sce){
  SGL_result_s <- foreach(r = 1:rep) %dopar% {
    if (s<=4) {SGL_manual_func(outcome_df = sim_outcome[[s]][[r]],
                               bm_df = sim_baseline[[s]][[r]])}
    else if (s==5) {SGL_manual_func(outcome_df = sim_outcome[[4]][[r]],
                                    bm_df = sim_baseline[[5]][[r]])}
  }
  SGL_result <- append(SGL_result, list(SGL_result_s))
}

# Run SGL for each scenario and each rep on data simulated with Model 2
SGL_result_quad <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has 5 rows and 3 cols for the CV-picked SGL beta coefficients
for (s in 1:num_sce){
  SGL_result_quad_s <- foreach(r = 1:rep) %dopar% {
    if (s<=4) {SGL_manual_func(outcome_df = sim_outcome_quad[[s]][[r]],
                               bm_df = sim_baseline_quad[[s]][[r]])}
    else if (s==5) {SGL_manual_func(outcome_df = sim_outcome_quad[[4]][[r]],
                                    bm_df = sim_baseline_quad[[5]][[r]])}
  }
  SGL_result_quad <- append(SGL_result_quad, list(SGL_result_quad_s))
}

# Run SGL for each scenario and each rep on data simulated with Model 3
SGL_result_splines <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has 5 rows and 3 cols for the CV-picked SGL beta coefficients
for (s in 1:num_sce){
  SGL_result_splines_s <- foreach(r = 1:rep) %dopar% {
    if (s<=4) {SGL_manual_func(outcome_df = sim_outcome_splines[[s]][[r]],
                               bm_df = sim_baseline_splines[[s]][[r]])}
    else if (s==5) {SGL_manual_func(outcome_df = sim_outcome_splines[[4]][[r]],
                                    bm_df = sim_baseline_splines[[5]][[r]])}
  }
  SGL_result_splines <- append(SGL_result_splines, list(SGL_result_splines_s))
}

# Clean results (mean and sd of beta coefficients, CV-picked lambda and lldiff) from SGL_manual_func()
# and count number of non-zero beta coefficients
SGL_result_clean <- clean_func(SGL_result)
View(SGL_result_clean$beta_clean)
View(SGL_result_clean$index_lambda)

SGL_result_quad_clean <- clean_func(SGL_result_quad)
View(SGL_result_quad_clean$beta_clean)
View(SGL_result_quad_clean$index_lambda)

SGL_result_splines_clean <- clean_func(SGL_result_splines)
View(SGL_result_splines_clean$beta_clean)
View(SGL_result_splines_clean$index_lambda)



