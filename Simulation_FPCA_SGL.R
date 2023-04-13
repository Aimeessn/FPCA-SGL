# Run our proposed FPCA-SGL method on simulated data (3 models and 5 scenarios as listed in Table 2)
# Calculate TPR and FPR

library(dplyr)
library(tidyverse)
library(fdapace)
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
  library(fdapace)
})
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()



################################################################################################################################
# Conduct FPCA for each biomarker

# Function to FPCA (with k=3) on each of the biomarker for each repetition
# Arguments: a dataframe "bm_df" for simulated observed biomarker trajectory
# Output: a list of `num_sce` lists (one for each simulation scenario),
#         each list is a list of `rep` dfs,
#         each df has `sample_size` rows and 13 cols (ID + 3FPCs for 4bms)
FPC3 <- function(bm_df){
  # Create a list of dfs where each df is a long format of one biomarker's measurements
  bm_group = data.frame(bm = c("bm1", "bm2", "bm3", "bm4"), group = c(1,1,2,2))
  bm_list <- list()
  for (i in 1:nrow(bm_group)){
    temp <- bm_df %>% dplyr::select(ID, Time, all_of(bm_group$bm[i]))
    bm_list[[i]] <- temp[complete.cases(temp),]
  }
  names(bm_list) <- bm_group$bm
  
  # Run FPCA on each biomarker
  fpca <- list()
  for (i in 1:length(bm_list)){
    bm_i <- bm_list[[i]]
    colnames(bm_i)[3] <- "bm"
    ID <- unique(bm_i$ID)
    yList <- list()
    tList <- list()
    for (j in 1:length(ID)){
      yList[[j]] <- bm_i$bm[bm_i$ID==j]
      tList[[j]] <- bm_i$Time[bm_i$ID==j]
    }
    fpca[[i]] <- FPCA(Ly=yList, Lt=tList, optns=list(methodSelectK = 3))
  }
  
  # Select k for the number of FPCs we will use for each biomarker
  cumFVE_1st <- c()
  cumFVE_2nd <- c()
  cumFVE_3rd <- c()
  selected_k <- c()
  for (i in 1:length(fpca)){
    cumFVE_1st[i] <- fpca[[i]]$cumFVE[1]
    cumFVE_2nd[i] <- fpca[[i]]$cumFVE[2]
    cumFVE_3rd[i] <- fpca[[i]]$cumFVE[3]
    selected_k[i] <- length(fpca[[i]]$cumFVE)
  }
  
  # Consolidate the first k=3 FPCs for each biomarker
  k = 3
  bm_fpc <- data.frame(ID = seq(1, sample_size))
  for (i in 1:length(fpca)){
    bm_fpc_i <- cbind(unique(bm_list[[i]]$ID), as.data.frame(fpca[[i]]$xiEst[,1:k]))
    colnames(bm_fpc_i) <- c("ID", paste(names(bm_list)[i], "_1st", sep=""),
                            paste(names(bm_list)[i], "_2nd", sep=""),
                            paste(names(bm_list)[i], "_3rd", sep=""))
    bm_fpc <- bm_fpc %>% left_join(bm_fpc_i, by = "ID")
  }
  return(bm_fpc)
}

# Run FPCA on data simulated with Model 1
sim_fpc <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has `sample_size` rows and 13 cols (ID + 3FPCs for 4bms)
for (s in 1:num_sce){
  sim_fpc_s <- foreach(r = 1:rep) %dopar% {
    FPC3(bm_df = sim_bm_obs[[s]][[r]])
  }
  sim_fpc <- append(sim_fpc, list(sim_fpc_s))
}

# Run FPCA on data simulated with Model 2
sim_fpc_quad <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has `sample_size` rows and 13 cols (ID + 3FPCs for 4bms)
for (s in 1:num_sce){
  sim_fpc_quad_s <- foreach(r = 1:rep) %dopar% {
    FPC3(bm_df = sim_bm_obs_quad[[s]][[r]])
  }
  sim_fpc_quad <- append(sim_fpc_quad, list(sim_fpc_quad_s))
}

# Run FPCA on data simulated with Model 3
sim_fpc_splines <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has `sample_size` rows and 13 cols (ID + 3FPCs for 4bms)
for (s in 1:num_sce){
  sim_fpc_splines_s <- foreach(r = 1:rep) %dopar% {
    FPC3(bm_df = sim_bm_obs_splines[[s]][[r]])
  }
  sim_fpc_splines <- append(sim_fpc_splines, list(sim_fpc_splines_s))
}



################################################################################################################################
# Conduct SGL on the first 3 FPCs, with pre-picked alpha=0.05 (~GL), and pick lambda using CV

# Main function to perform SGL with pre-picked alpha=0.05 (~GL), and pick lambda using CV
# The range of 50 potential lambda values is defined as:
#   upper bound being the smallest lambda that shrink everything to 0, lower bound
#   lower bound is upper bound / 100
# lambda_min: the lambda that gives the minimum CV lldiff
# lambda_1se: the largest lambda (more shrinkage) that gives the CV lldiff within 1 standard error of the minimum
# Arguments: a dataframe "outcome_df" for outcome data simulated under different scenario
#            a dataframe "bm_df" for simulated biomarkers data (first 3 FPCs)
# Output: a list of `num_sce` lists (one for each simulation scenario),
#         each list is a list of `rep` dfs,
#         each df has p+2 (p is the number of covariates, not including intercept, the other 2 is one for
#         lambda index and one for lambda) rows and 3 columns (variable, bm_group, beta_best)
SGL_manual_func <- function(outcome_df, bm_df){
  # Combine biomarker, demographics, and outcome data
  df <- outcome_df %>% dplyr::select(ID, time, status) %>%
    left_join(bm_df, by = "ID")
  
  # Standardize biomarker data
  df[c(4:15)] <- scale(df[c(4:15)])
  
  # Prepare for SGL
  SGL_x <- model.matrix(terms( ~ bm1_1st + bm1_2nd + bm1_3rd + bm2_1st + bm2_2nd + bm2_3rd +
                                 bm3_1st + bm3_2nd + bm3_3rd + bm4_1st + bm4_2nd + bm4_3rd, keep.order = TRUE),
                        df %>% dplyr::select(-ID, -time, -status))
  SGL_list <- list(x = SGL_x[,-1], time = df$time, status = df$status)
  SGL_group <- c(rep(1,6),rep(2,6))
  
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
  if (SGL_index != length(SGL_index)){SGL_index = lambda_1se_index} 
  
  # Clean final model and save beta coeffcients and CV-picked lambda and alpha
  variable = c(colnames(df)[4:15], "index", "lambda")
  bm_group = c(rep("Group 1",6), rep("Group 2",6), "index", "lambda")
  result <- data.frame(beta_best = beta_all[,SGL_index]) %>%
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
      beta[,r] <- result_list[[s]][[r]]$beta_best[-c(13:14)]
      if (result_list[[s]][[r]]$beta_best[13] == index) {index_numlargest = index_numlargest + 1}
      lambda_s[r] = result_list[[s]][[r]]$beta_best[14]
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
SGL_result <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has 12 rows and 3 cols for the CV-picked SGL beta coefficients
for (s in 1:num_sce){
  SGL_result_s <- foreach(r = 1:rep) %dopar% {
    if (s<=4) {SGL_manual_func(outcome_df = sim_outcome[[s]][[r]],
                               bm_df = sim_fpc[[s]][[r]])}
    else if (s==5) {SGL_manual_func(outcome_df = sim_outcome[[4]][[r]],
                                    bm_df = sim_fpc[[5]][[r]])}
  }
  SGL_result <- append(SGL_result, list(SGL_result_s))
}

# Run SGL for each scenario and each rep on data simulated with Model 2
SGL_result_quad <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has 12 rows and 3 cols for the CV-picked SGL beta coefficients
for (s in 1:num_sce){
  SGL_result_quad_s <- foreach(r = 1:rep) %dopar% {
    if (s<=4) {SGL_manual_func(outcome_df = sim_outcome_quad[[s]][[r]],
                               bm_df = sim_fpc_quad[[s]][[r]])}
    else if (s==5) {SGL_manual_func(outcome_df = sim_outcome_quad[[4]][[r]],
                                    bm_df = sim_fpc_quad[[5]][[r]])}
  }
  SGL_result_quad <- append(SGL_result_quad, list(SGL_result_quad_s))
}

# Run SGL for each scenario and each rep on data simulated with Model 3
SGL_result_splines <- list() #a list of `num_sce` lists (one for each simulation scenario), each list is a list of `rep` dfs, each df has 12 rows and 3 cols for the CV-picked SGL beta coefficients
for (s in 1:num_sce){
  SGL_result_splines_s <- foreach(r = 1:rep) %dopar% {
    if (s<=4) {SGL_manual_func(outcome_df = sim_outcome_splines[[s]][[r]],
                               bm_df = sim_fpc_splines[[s]][[r]])}
    else if (s==5) {SGL_manual_func(outcome_df = sim_outcome_splines[[4]][[r]],
                                    bm_df = sim_fpc_splines[[5]][[r]])}
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


