#### Description: Simple probabilistic state-transition model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array
library(parallel)
library(future.apply)
library(microbenchmark)  # to track time

# C++ code requires Rtools 
library(Rcpp)
library(RcppArmadillo)

# Clear workspace
rm(list = ls()) # clear objects from workspace

# Load functions from local scripts
source("f_model.R")
source("f_model_df_conversion.R")
source("f_model_matrix.R")
source("f_model_test.R")
source("f_model_test_rcpp.R")

# load loop function in C++
sourceCpp("f_propagate_states.cpp")

#### General ####
v_states <- c("Gezond", "Ziek", "Dood") #  vector of model health states
n_states <- length(v_states) # number of health states 
v_treatments <- c("Current_practice", "New_treatment") # vector of strategy names
n_treatments <- length(v_treatments)  # number of treatments
n_t <- 100 # model time horizon 
n_sim <- 5000 #00 # number of Monte-carlo simulations for probabilistic analyses
set.seed(12345) # set seed

#### Model inputs #### 
df_input <- data.frame( # open input parameter dataframe
  # transition probabilities 
  p_gezond_ziek = rbeta(n = n_sim, shape1 = 10, shape2 = 20), 
  p_gezond_dood = rbeta(n = n_sim, shape1 = 1, shape2 = 20), 
  p_ziek_gezond = rbeta(n = n_sim, shape1 = 5, shape2 = 20), 
  p_ziek_dood = rbeta(n = n_sim, shape1 = 10, shape2 = 20),
  # relative risk (difference between treatments)
  rr_gezond_ziek_t2_t1 = rlnorm(n = n_sim, meanlog = log(0.8), sdlog = 0.2),
  # health state utility (quality of life)
  u_gezond = rbeta(n = n_sim, shape1 = 75, shape2 = 100), 
  u_ziek = rbeta(n = n_sim, shape1 = 60, shape2 = 100), 
  u_dood = rep(x = 0, times = n_sim),
  # health state costs
  c_gezond = rgamma(n = n_sim, shape = 75, scale = 100), 
  c_ziek = rgamma(n = n_sim, shape = 100, scale = 100), 
  c_dood = rep(x = 0, times = n_sim)
) # close input parameter dataframe

#### Create matrix to store results #### 
m_out_1a <- m_out_2a <- m_out_3a <- m_out_4a <- m_out_1b <- m_out_2b <- m_out_3b <- m_out_4b <- m_out_1c <- m_out_2c <- m_out_3c <- m_out_4c <- m_out_4d <- m_out_4e <- matrix( 
  data = NA, 
  nrow = n_sim, 
  ncol = 4,
  dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments))) 
) # end matrix

#### Evaluate model calculation approaches with microbenchmark() #### 
microbenchmark(
  approach1 = {
    # 1a Run model calculations with for loop and df conversion for mapply
    for(x in 1:n_sim){ # loop to run model for each row of PSA inputs
      m_out_1a[x, ] <- f_model_df_conversion(x)
    } # close for loop
  },
  
  approach2 = {
    # 1b Run model calculations with for loop 
    for(x in 1:n_sim){ # loop to run model for each row of PSA inputs
      m_out_1b[x, ] <- f_model(x)
    } # close for loop
  },
  
  approach3 = {
    # 1c Run model calculations with for loop and matrix as input
    for(x in 1:n_sim){ # loop to run model for each row of PSA inputs
      m_out_1c[x,] <- f_model_matrix(as.matrix(df_input[x,]))
    } # close for loop
  },
  
  approach4 = {
    # 2a Run model calculations with apply and df conversion for mapply
    m_out_2a <- t(apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_df_conversion(x)
    }))
  },
  
  approach5 = {
    # 2b Run model calculations with apply
    m_out_2b <- t(apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model(x)
    }))
  },
  
  approach6 = {
    # 2c Run model calculations with apply and matrix as input
    m_out_2c <- t(apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_matrix(as.matrix(df_input[x,]))
    }))
  },
  
  approach7 = {
    # 3a Run model with parallel processing and df conversion for mapply
    cl <- makeCluster(detectCores()) # Create a cluster with the number of cores available in your machine
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input")) # Export necessary variables to the cluster
    m_out_3a <- parLapply(cl, 1:n_sim, f_model_df_conversion)
    m_out_3a <- do.call(rbind, m_out_3a) # combine list into a matrix
    stopCluster(cl)
  },
  
  approach8 = {
    # 3b Run model with parallel processing 
    cl <- makeCluster(detectCores()) # Create a cluster with the number of cores available in your machine
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input")) # Export necessary variables to the cluster
    m_out_3b <- parLapply(cl, 1:n_sim, f_model)
    m_out_3b <- do.call(rbind, m_out_3b) # combine list into a matrix
    stopCluster(cl)
  },
  
  approach9 = {
    # 3c Run model with parallel processing and matrix as input
    cl <- makeCluster(detectCores()) # Create a cluster with the number of cores available in your machine
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_matrix")) # Export necessary variables to the cluster
    m_out_3c <- parLapply(cl, 1:n_sim, function(x) f_model_matrix(as.numeric(df_input[x, ])))
    m_out_3c <- do.call(rbind, m_out_3c) # combine list into a matrix
    stopCluster(cl)
  },
  
  approach10 = {
    # 4a Run model with calculations with future apply and df conversion for mapply
    plan(multisession) ## Run in parallel on local computer
    m_out_4a <- t(future_apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_df_conversion(x)
    }))
    plan(sequential) # Stop running in parallel on local computer
  },
  
  approach11 = {
    # 4b Run model with calculations with future apply
    plan(multisession) ## Run in parallel on local computer
    m_out_4b <- t(future_apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model(x)
    }))
    plan(sequential) # Stop running in parallel on local computer
  },
  
  approach12 = {
    # 4c Run model with calculations with future apply + inputs in matrix
    plan(multisession) ## Run in parallel on local computer
    m_out_4c <- t(future_apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_matrix(as.matrix(df_input[x,]))
    }))
    plan(sequential) # Stop running in parallel on local computer
  },
  
  approach13 = {
    # 4 Run model with calculations with future apply + test
    plan(multisession) ## Run in parallel on local computer
    m_out_4d <- future_sapply(1:n_sim, function(x) {
      f_model_test(as.numeric(df_input[x, ]))
    })
    plan(sequential) # Stop running in parallel on local computer
  },
  
  approach14 = {
    # 4 Run model with calculations with future apply + test + using rccp
    plan(multisession) ## Run in parallel on local computer
    m_out_4e <- future_sapply(1:n_sim, function(x) {
      f_model_test_rcpp(as.numeric(df_input[x, ]))
    })
    plan(sequential) # Stop running in parallel on local computer
  },
  
  times = 50 # how often each approach should be evaluated
)

summary(m_out_1a == m_out_2a) 
summary(m_out_2a == m_out_3a) 
summary(m_out_3a == m_out_4a)
summary(m_out_4a == m_out_1b) 
summary(m_out_1b == m_out_2b)
summary(m_out_2b == m_out_3b)
summary(m_out_3b == m_out_4b)
summary(m_out_4b == m_out_1c)
summary(m_out_1c == m_out_2c)
summary(m_out_2c == m_out_3c)
summary(m_out_3c == m_out_4c)
summary(m_out_4c == t(m_out_4d))
