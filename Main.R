#### Description: Simple Probabilistic State-Transition Model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array ####

# General settings
options(scipen = 999, max.print = 10000, digits = 4)

# Load and install necessary libraries
required_packages <- c(
  "parallel", "future.apply", "microbenchmark", "Rcpp", "RcppArmadillo", "docstring"
)
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
suppressPackageStartupMessages(lapply(required_packages, require, character.only = TRUE))

# Clear workspace
rm(list = ls())

# Load custom functions
source("f_model_a.R")
source("f_model_b.R")
source("f_model_c.R")
source("f_model_d.R")
source("f_model_e.R")

# Load Rcpp function for state propagation
sourceCpp("f_propagate_states.cpp")

#### General Setup ####
v_states <- c("Healthy", "Sick", "Death")               # Vector of model health states
n_states <- length(v_states)                            # Number of health states
v_treatments <- c("Current_practice", "New_treatment")  # Strategy names
n_treatments <- length(v_treatments)                    # Number of treatments
n_t <- 100                                              # Model time horizon
n_sim <- 5000                                           # Number of Monte Carlo simulations
n_iterations <- 20                                      # Define the number of iterations for benchmark
set.seed(12345)                                         # Set random seed

#### Model Inputs ####
# Create a dataframe for probabilistic sensitivity analysis (PSA) inputs
df_input <- data.frame(
  # Transition probabilities
  p_healthy_sick = rbeta(n = n_sim, shape1 = 10, shape2 = 20),
  p_healthy_death = rbeta(n = n_sim, shape1 = 1, shape2 = 20),
  p_sick_healthy = rbeta(n = n_sim, shape1 = 5, shape2 = 20),
  p_sick_death = rbeta(n = n_sim, shape1 = 10, shape2 = 20),
  # Relative risk for new treatment
  rr_healthy_sick_t2_t1 = rlnorm(n = n_sim, meanlog = log(0.8), sdlog = 0.2),
  # Utilities
  u_healthy = rbeta(n = n_sim, shape1 = 75, shape2 = 100),
  u_sick = rbeta(n = n_sim, shape1 = 60, shape2 = 100),
  u_death = rep(0, n_sim),
  # Costs
  c_healthy = rgamma(n = n_sim, shape = 75, scale = 100),
  c_sick = rgamma(n = n_sim, shape = 100, scale = 100),
  c_death = rep(0, n_sim)
)

# f_model_a(df_input[1, ])
# f_model_b(df_input[1, ])
# f_model_c(as.matrix(df_input[1, ]))
# f_model_d(df_input[1, ])
# f_model_e(df_input[1, ])

#### Create Matrices for Results ####
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

#### Evaluate Model Calculation Approaches ####
df_benchmark_results <- microbenchmark(
  approach1 = {  # Sequential loop: f_model_a
    for (x in 1:n_sim) l_result_matrices[[1]][x, ] <- f_model_a(df_input[x, ])
  },
  
  approach2 = {  # Sequential loop: f_model_b
    for (x in 1:n_sim) l_result_matrices[[2]][x, ] <- f_model_b(df_input[x, ])
  },
  
  approach3 = {  # Sequential loop: f_model_c
    for (x in 1:n_sim) l_result_matrices[[3]][x, ] <- f_model_c(as.matrix(df_input[x, ]))
  },
  
  approach4 = {  # Sequential loop: f_model_d
    for (x in 1:n_sim) l_result_matrices[[4]][x, ] <- f_model_d(df_input[x, ])
  },
  
  approach5 = {  # Sequential loop: f_model_e
    for (x in 1:n_sim) l_result_matrices[[5]][x, ] <- f_model_e(df_input[x, ])
  },
  
  approach6 = {  # Apply-based approach: f_model_a
    l_result_matrices[[6]] <- t(apply(matrix(seq_len(nrow(df_input))), 1, function(x) f_model_a(df_input[x, ])))
  },
  
  approach7 = {  # Apply-based approach: f_model_b
    l_result_matrices[[7]] <- t(apply(matrix(seq_len(nrow(df_input))), 1, function(x) f_model_b(df_input[x, ])))
  },
  
  approach8 = {  # Apply-based approach: f_model_c
    l_result_matrices[[8]] <- t(apply(matrix(seq_len(nrow(df_input))), 1, function(x) f_model_c(as.matrix(df_input[x, ]))))
  },
  
  approach9 = {  # Apply-based approach: f_model_d
    l_result_matrices[[9]] <- t(apply(matrix(seq_len(nrow(df_input))), 1, function(x) f_model_d(df_input[x, ])))
  },
  
  approach10 = {  # Apply-based approach: f_model_e
    l_result_matrices[[10]] <- t(apply(matrix(seq_len(nrow(df_input))), 1, function(x) f_model_e(df_input[x, ])))
  },
  
  approach11 = {  # Parallel: f_model_a
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_a"))
    l_result_matrices[[11]] <- do.call(rbind, parLapply(cl, 1:n_sim, function(x) f_model_a(df_input[x, ])))
    stopCluster(cl)
  },
  
  approach12 = {  # Parallel: f_model_b
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_b"))
    l_result_matrices[[12]] <- do.call(rbind, parLapply(cl, 1:n_sim, function(x) f_model_b(df_input[x, ])))
    stopCluster(cl)
  },
  
  approach13 = {  # Parallel: f_model_c
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_c"))
    l_result_matrices[[13]] <- do.call(rbind, parLapply(cl, 1:n_sim, function(x) f_model_c(as.matrix(df_input[x, ]))))
    stopCluster(cl)
  },
  
  approach14 = {  # Parallel: f_model_d
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_d"))
    l_result_matrices[[14]] <- do.call(rbind, parLapply(cl, 1:n_sim, function(x) f_model_d(df_input[x, ])))
    stopCluster(cl)
  },
  
  approach15 = {  # Parallel: f_model_e
    cl <- makeCluster(detectCores())
    clusterEvalQ(cl, Rcpp::sourceCpp("f_propagate_states.cpp"))
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_e"))
    l_result_matrices[[15]] <- do.call(rbind, parLapply(cl, 1:n_sim, function(x) f_model_e(df_input[x, ])))
    stopCluster(cl)
  },
  
  approach16 = {  # Future multisession: f_model_a
    plan(multisession)
    l_result_matrices[[16]] <- t(future_sapply(1:n_sim, function(x) f_model_a(df_input[x, ]), future.seed = TRUE))
    plan(sequential)
  },
  
  approach17 = {  # Future multisession: f_model_b
    plan(multisession)
    l_result_matrices[[17]] <- t(future_sapply(1:n_sim, function(x) f_model_b(df_input[x, ]), future.seed = TRUE))
    plan(sequential)
  },
  
  approach18 = {  # Future multisession: f_model_c
    plan(multisession)
    l_result_matrices[[18]] <- t(future_sapply(1:n_sim, function(x) f_model_c(as.matrix(df_input[x, ])), future.seed = TRUE))
    plan(sequential)
  },
  
  approach19 = {  # Future multisession: f_model_d
    plan(multisession)
    l_result_matrices[[19]] <- t(future_sapply(1:n_sim, function(x) f_model_d(df_input[x, ]), future.seed = TRUE))
    plan(sequential)
  },
  
  approach20 = {  # Future multisession: f_model_e
    plan(multisession)
    l_result_matrices[[20]] <- t(future_sapply(1:n_sim, function(x) {
      Rcpp::sourceCpp("f_propagate_states.cpp")
      f_model_e(df_input[x, ])
    }, future.seed = TRUE))
    plan(sequential)
  },
  times = n_iterations
)

df_summary <- as.data.frame(summary(df_benchmark_results))
m_results <- matrix(data = df_summary$mean, nrow = 5, ncol = 4, byrow = FALSE,
                  dimnames = list(c("model a", "model b", "model c", "model d", "model e"), 
                                  c("Sequential loop", "Apply-based approach", "Parallel", "Future multisession")))

m_results
colSums(m_results)
rowSums(m_results)
m_results[ , 1] - m_results[ , 2]
m_results[ , 3] - m_results[ , 4]

write.csv(df_summary, file = paste0("benchmark n_sim = ", n_sim, " n_iterations = ", n_iterations, ".csv"), row.names = FALSE)
write.csv(m_results, file = paste0("benchmark mean n_sim = ", n_sim, " n_iterations = ", n_iterations, ".csv"), row.names = FALSE)

# sum(df_summary$mean) / 60 # time (minutes) per iteration

#### Compare Result Matrices ####
# Pairwise comparison of matrices in the result list
l_comparison_results <- lapply(seq_along(l_result_matrices), function(i) {
  if (i == length(l_result_matrices)) {
    l_result_matrices[[i]] - l_result_matrices[[1]]  # Compare last with first
  } else {
    l_result_matrices[[i]] - l_result_matrices[[i + 1]]  # Compare with next
  }
})

# Apply summary to all matrices in l_comparison_results
l_summaries <- lapply(l_comparison_results, summary)

# View all summaries at once (prints each one sequentially)
l_summaries
