#### Description: Simple Probabilistic State-Transition Model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array, l = list ####

# Load packages ----------------------------------------------------------------

# Check if pacman package is installed. If not, install it

if(!require(pacman)) install.packages("pacman")

# List needed packages

v_packages <- c("bench", # for benchmarking
                "parallel", # for parallel processing
                "future.apply", # for parallel processing
                "Rcpp", # for C++ code
                "RcppArmadillo", # for C++ code
                "docstring", # for package-like documentation
                "tidyverse") # for data manipulation

# Load packages

pacman::p_load(v_packages, character.only = TRUE)


# Settings ----------------------------------------------------------------

# General settings
options(scipen = 999, max.print = 10000)

# Define the number of iterations per benchmark
n_iterations <- 2

# Load custom functions
source("f_model_a.R") # check function using: docstring(f_model_a)
source("f_model_b.R") # check function using: docstring(f_model_b)
source("f_model_c.R") # check function using: docstring(f_model_c)
source("f_model_d.R") # check function using: docstring(f_model_d)
source("f_model_e.R") # check function using: docstring(f_model_e)

# Load Rcpp function for state propagation
sourceCpp("f_propagate_states.cpp")


#### General model set-up ####
v_states <- c("Healthy", "Sick", "Death") #  vector of model health states
n_states <- length(v_states) # number of health states 
v_treatments <- c("Current_practice", "New_treatment") # vector of strategy names
n_treatments <- length(v_treatments)  # number of treatments
n_t <- 100 # model time horizon 
n_sim <- 50  # Number of Monte Carlo simulations
set.seed(134340) # set seed

#### Model inputs #### 
df_input <- data.frame( # open input parameter dataframe
  # transition probabilities 
  p_healthy_sick = rbeta(n = n_sim, shape1 = 10, shape2 = 20), 
  p_healthy_death = rbeta(n = n_sim, shape1 = 1, shape2 = 20), 
  p_sick_healthy = rbeta(n = n_sim, shape1 = 5, shape2 = 20), 
  p_sick_death = rbeta(n = n_sim, shape1 = 10, shape2 = 20),
  # relative risk (difference between treatments)
  rr_healthy_sick_t2_t1 = rlnorm(n = n_sim, meanlog = log(0.8), sdlog = 0.2),
  # health state utility (quality of life)
  u_healthy = rbeta(n = n_sim, shape1 = 75, shape2 = 100), 
  u_sick = rbeta(n = n_sim, shape1 = 60, shape2 = 100), 
  u_death = rep(x = 0, times = n_sim),
  # health state costs
  c_healthy = rgamma(n = n_sim, shape = 75, scale = 100), 
  c_sick = rgamma(n = n_sim, shape = 100, scale = 100), 
  c_death = rep(x = 0, times = n_sim)
) # close input parameter dataframe

#### Create matrix to store results #### 
# Use a list to store result matrices
l_result_matrices <- setNames(vector("list", 20), paste0("m_out_", 1:20))

# Initialize each matrix in the list
for (i in seq_along(l_result_matrices)) {
  l_result_matrices[[i]] <- matrix(
    data = NA,
    nrow = n_sim,
    ncol = 2 * n_treatments,
    dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments)))
  )
}

# Benchmark each approach
benchmark_results <- bench::mark(
  # loop approach
  a_loop = 
    # Sequential loop: f_model_a
    for (x in 1:n_sim) l_result_matrices[[1]][x, ] <- f_model_a(df_input[x, ]),
  b_loop = 
    # Sequential loop: f_model_b
    for (x in 1:n_sim) l_result_matrices[[2]][x, ] <- f_model_b(df_input[x, ]),
  c_loop =
     # Sequential loop: f_model_c
     for (x in 1:n_sim) l_result_matrices[[3]][x, ] <- f_model_c(df_input[x, ]),
  d_loop =
     # Sequential loop: f_model_d
     for (x in 1:n_sim) l_result_matrices[[4]][x, ] <- f_model_d(df_input[x, ]),
  e_loop =
    # Sequential loop: f_model_e
    for (x in 1:n_sim) l_result_matrices[[5]][x, ] <- f_model_e(df_input[x, ]),
  # Apply-based approach
  # Apply-based: f_model_a
  a_apply = 
    l_result_matrices[[6]] <- t(apply(matrix(seq_len(nrow(df_input))),
                                      1, function(x) f_model_a(df_input[x, ]))),
  # Apply-based: f_model_b
  b_apply = 
    l_result_matrices[[7]] <- t(apply(matrix(seq_len(nrow(df_input))),
                                      1, function(x) f_model_b(df_input[x, ]))),
  # Apply-based: f_model_c
  c_apply =
    l_result_matrices[[8]] <- t(apply(matrix(seq_len(nrow(df_input))),
                                      1, function(x) f_model_c(df_input[x, ]))),
  # Apply-based: f_model_d
  d_apply =
    l_result_matrices[[9]] <- t(apply(matrix(seq_len(nrow(df_input))),
                                      1, function(x) f_model_d(df_input[x, ]))),
  # Apply-based: f_model_e
  e_apply =
    l_result_matrices[[10]] <- t(apply(matrix(seq_len(nrow(df_input))),
                                      1, function(x) f_model_e(df_input[x, ]))),
  # Parallel-based approach
  # Parallel: f_model_a
  a_parallel = {
    # Parallel-based: f_model_a
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments",
                        "n_t", "df_input", "f_model_a"))
    l_result_matrices[[11]] <- do.call(rbind, parLapply(cl,
                                                        1:n_sim,
                                                        function(x) f_model_a(
                                                          df_input[x, ])))
    stopCluster(cl)},
  # Parallel: f_model_b
  b_parallel = {
    # Parallel-based: f_model_b
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments",
                        "n_t", "df_input", "f_model_b"))
    l_result_matrices[[12]] <- do.call(rbind, parLapply(cl,
                                                        1:n_sim,
                                                        function(x) f_model_b(
                                                          df_input[x, ])))
    stopCluster(cl)},
  # Parallel: f_model_c
  c_parallel = {
    # Parallel-based: f_model_c
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments",
                        "n_t", "df_input", "f_model_c"))
    l_result_matrices[[13]] <- do.call(rbind, parLapply(cl,
                                                        1:n_sim,
                                                        function(x) f_model_c(
                                                          df_input[x, ])))
    stopCluster(cl)},
  # Parallel: f_model_d
  d_parallel = {
    # Parallel-based: f_model_d
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments",
                        "n_t", "df_input", "f_model_d"))
    l_result_matrices[[14]] <- do.call(rbind, parLapply(cl,
                                                        1:n_sim,
                                                        function(x) f_model_d(
                                                          df_input[x, ])))
    stopCluster(cl)},
  # Parallel: f_model_e
  e_parallel = {
    # Parallel-based: f_model_e
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments",
                        "n_t", "df_input", "f_model_e"))
    l_result_matrices[[15]] <- do.call(rbind, parLapply(cl,
                                                        1:n_sim,
                                                        function(x) f_model_e(
                                                          df_input[x, ])))
    stopCluster(cl)},
  
  iterations = n_iterations,
  check = F
)

benchmark_results

# qualitative scores
# 
# v_qual_scores_rater1 <- c(5, 3)
# v_qual_scores_rater2 <- c(4, 3)
# 
# v_qual_score_mean <- (v_qual_scores_rater1 + v_qual_scores_rater2) / 2
# 
# benchmark_results$qual_mean <- v_qual_score_mean
# 
# # weights
# w1 <- 0.3
# w2 <- 0.3
# w3 <- 1-sum(w1, w2)
# 
# # overall score
# 
# benchmark_results$overall_score <- w1 * as.numeric(benchmark_results$total_time) +
#   w2 * as.numeric(benchmark_results$mem_alloc) + w3 * as.numeric(benchmark_results$qual_mean)
# 
# # Load functions from local scripts
# source("f_model.R")
# source("f_model_df_conversion.R")
# source("f_model_matrix.R")
