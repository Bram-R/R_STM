#### Description: Simple probabilistic state-transition model ####
#### Conventions: n = single number, v = vector, df = dataframe, m = matrix, a = array
options(scipen = 999) # setting for scientific notation
options(max.print = 10000) # setting for maximum output to display
options(digits = 4) # setting number of digits to display

library(parallel)
library(future.apply)
library(microbenchmark)  # to track time

# C++ code requires Rtools 
library(Rcpp)
library(RcppArmadillo)

# Clear workspace
rm(list = ls()) # clear objects from workspace

# Load functions from local scripts
source("f_model_a.R")
source("f_model_b.R")
source("f_model_c.R")
source("f_model_d.R")
source("f_model_e.R")

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

# validation
# f_model_a(df_input[1,])
# f_model_b(df_input[1,])
# f_model_c(as.matrix(df_input[1,]))
# f_model_d(df_input[1,])
# f_model_e(df_input[1,])

#### Create matrix to store results #### 
m_out_1a <- m_out_2a <- m_out_3a <- m_out_4a <- 
  m_out_1b <- m_out_2b <- m_out_3b <- m_out_4b <- 
  m_out_1c <- m_out_2c <- m_out_3c <- m_out_4c <- 
  m_out_1d <- m_out_2d <- m_out_3d <- m_out_4d <- 
  m_out_1e <- m_out_2e <- m_out_3e <- m_out_4e <- 
  matrix( 
    data = NA, 
    nrow = n_sim, 
    ncol = 4,
    dimnames = list(1:n_sim, c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments))) 
  ) # end matrix

#### Evaluate model calculation approaches with microbenchmark() #### 
microbenchmark(
  approach1 = {
    # 1a 
    for(x in 1:n_sim){ # loop to run model for each row of PSA inputs
      m_out_1a[x, ] <- f_model_a(df_input[x,])
    } # close for loop
  },
  
  approach2 = {
    # 1b 
    for(x in 1:n_sim){ # loop to run model for each row of PSA inputs
      m_out_1b[x, ] <- f_model_b(df_input[x,])
    } # close for loop
  },
  
  approach3 = {
    # 1c 
    for(x in 1:n_sim){ # loop to run model for each row of PSA inputs
      m_out_1c[x,] <- f_model_c(as.matrix(df_input[x,]))
    } # close for loop
  },
  
  approach4 = {
    # 1d 
    for(x in 1:n_sim){ # loop to run model for each row of PSA inputs
      m_out_1d[x, ] <- f_model_d(df_input[x,])
    } # close for loop
  },
  
  approach5 = {
    # 1e 
    for(x in 1:n_sim){ # loop to run model for each row of PSA inputs
      m_out_1e[x, ] <- f_model_e(df_input[x,])
    } # close for loop
  },
  
  approach6 = {
    # 2a 
    m_out_2a <- t(apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_a(df_input[x,])
    }))
  },
  
  approach7 = {
    # 2b 
    m_out_2b <- t(apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_b(df_input[x,])
    }))
  },
  
  approach8 = {
    # 2c 
    m_out_2c <- t(apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_c(as.matrix(df_input[x,]))
    }))
  },
  
  approach9 = {
    # 2d
    m_out_2d <- t(apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_d(df_input[x,])
    }))
  },
  
  approach10 = {
    # 2e
    m_out_2e <- t(apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_e(df_input[x,])
    }))
  },
  
  approach11 = {
    # 3a 
    cl <- makeCluster(detectCores()) # Create a cluster with the number of cores available in your machine
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_a")) # Export necessary variables to the cluster
    m_out_3a <- parLapply(cl, 1:n_sim, function(x) f_model_a(df_input[x, ]))
    m_out_3a <- do.call(rbind, m_out_3a) # combine list into a matrix
    stopCluster(cl)
  },
  
  approach12 = {
    # 3b 
    cl <- makeCluster(detectCores()) # Create a cluster with the number of cores available in your machine
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_b")) # Export necessary variables to the cluster
    m_out_3b <- parLapply(cl, 1:n_sim, function(x) f_model_b(df_input[x, ]))
    m_out_3b <- do.call(rbind, m_out_3b) # combine list into a matrix
    stopCluster(cl)
  },
  
  approach13 = {
    # 3c 
    cl <- makeCluster(detectCores()) # Create a cluster with the number of cores available in your machine
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_c")) # Export necessary variables to the cluster
    m_out_3c <- parLapply(cl, 1:n_sim, function(x) f_model_c(as.numeric(df_input[x, ])))
    m_out_3c <- do.call(rbind, m_out_3c) # combine list into a matrix
    stopCluster(cl)
  },
  
  approach14 = {
    # 3d 
    cl <- makeCluster(detectCores()) # Create a cluster with the number of cores available in your machine
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_d")) # Export necessary variables to the cluster
    m_out_3d <- parLapply(cl, 1:n_sim, function(x) f_model_d(df_input[x, ]))
    m_out_3d <- do.call(rbind, m_out_3d) # combine list into a matrix
    stopCluster(cl)
  },
  
  approach15 = {
    # 3e
    cl <- makeCluster(detectCores()) # Create a cluster with the number of cores available in your machine
    clusterEvalQ(cl, {Rcpp::sourceCpp("f_propagate_states.cpp")})  # Rcpp functions compiled using sourceCpp() are available in the main R session but are not automatically exported to parallel workers, which typically start as separate R processes. Hence it requires explicit loading in workers (or compiling a package)
    clusterExport(cl, c("v_states", "n_states", "v_treatments", "n_treatments", "n_t", "df_input", "f_model_e")) # Export necessary variables to the cluster
    m_out_3e <- parLapply(cl, 1:n_sim, function(x) f_model_e(df_input[x, ]))
    m_out_3e <- do.call(rbind, m_out_3e) # combine list into a matrix
    stopCluster(cl)
  },
  
  approach16 = {
    # 4a 
    plan(multisession) ## Run in parallel on local computer
    m_out_4a <- t(future_apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_a(df_input[x,])
    }))
    plan(sequential) # Stop running in parallel on local computer
  },
  
  approach17 = {
    # 4b 
    plan(multisession) ## Run in parallel on local computer
    m_out_4b <- t(future_apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_b(df_input[x,])
    }))
    plan(sequential) # Stop running in parallel on local computer
  },
  
  approach18 = {
    # 4c 
    plan(multisession) ## Run in parallel on local computer
    m_out_4c <- t(future_apply(X = matrix(seq_len(nrow(df_input))), MARGIN = 1, function(x) {
      f_model_c(as.matrix(df_input[x,]))
    }))
    plan(sequential) # Stop running in parallel on local computer
  },
  
  approach19 = {
    # 4d 
    plan(multisession) ## Run in parallel on local computer
    m_out_4d <- t(future_sapply(1:n_sim, function(x) {
      f_model_d(df_input[x, ])
    }))
    plan(sequential) # Stop running in parallel on local computer
  },
  
  approach20 = {
    # 4e
    plan(multisession) ## Run in parallel on local computer
    m_out_4e <- t(future_sapply(1:n_sim, function(x) {
      Rcpp::sourceCpp("f_propagate_states.cpp") # Rcpp functions compiled using sourceCpp() are available in the main R session but are not automatically exported to parallel workers, which typically start as separate R processes. Hence it requires explicit loading in workers (or compiling a package) 
      f_model_e(df_input[x, ])
    }, future.seed = TRUE))  # Ensures reproducible and safe random number generation
    plan(sequential) # Stop running in parallel on local computer
  },
  
  times = 5#0 # how often each approach should be evaluated
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
summary(m_out_4c == m_out_1d)
summary(m_out_1d == m_out_2d)
summary(m_out_2d == m_out_3d)
summary(m_out_3d == t(m_out_4d))
summary(t(m_out_4d) == m_out_1e)
summary(m_out_1e == m_out_2e)
summary(m_out_2e == m_out_3e)
summary(m_out_3e == t(m_out_4e))
summary(t(m_out_4e) == m_out_1a)
