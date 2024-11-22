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
source("f_model_a.R")
source("f_model_b.R")
source("f_model_c.R")
source("f_model_d.R")
source("f_model_e.R")

# load loop function in C++
sourceCpp("f_propagate_states.cpp")

#### General ####
v_states <- c("Healthy", "Sick", "Death") #  vector of model health states
n_states <- length(v_states) # number of health states 
v_treatments <- c("Current_practice", "New_treatment") # vector of strategy names
n_treatments <- length(v_treatments)  # number of treatments
n_t <- 100 # model time horizon 
n_sim <- 5 #0#00 #00 # number of Monte-carlo simulations for probabilistic analyses
set.seed(12345) # set seed

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


f_model_a(df_input[1,])
f_model_b(df_input[1,])
f_model_c(as.matrix(df_input[1,]))
f_model_d(df_input[1,])
f_model_e(df_input[1,])

