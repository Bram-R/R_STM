#### Model calculations by ChatGPT ####

f_model_test <- function(params) {
  # Extract parameters for clarity
  p_gezond_ziek <- params[1]
  p_gezond_dood <- params[2]
  p_ziek_gezond <- params[3]
  p_ziek_dood <- params[4]
  rr_gezond_ziek_t2_t1 <- params[5]
  u_gezond <- params[6]
  u_ziek <- params[7]
  u_dood <- params[8]
  c_gezond <- params[9]
  c_ziek <- params[10]
  c_dood <- params[11]
  
  # Transition matrix setup for Current Practice
  P_cp <- matrix(0, nrow = n_states, ncol = n_states)
  P_cp[1, ] <- c(1 - p_gezond_ziek - p_gezond_dood, p_gezond_ziek, p_gezond_dood)  # Gezond
  P_cp[2, ] <- c(p_ziek_gezond, 1 - p_ziek_gezond - p_ziek_dood, p_ziek_dood)      # Ziek
  P_cp[3, 3] <- 1  # Dood
  
  # Transition matrix setup for New Treatment
  P_nt <- P_cp
  P_nt[1, 2] <- p_gezond_ziek * rr_gezond_ziek_t2_t1  # Adjust transition rate to Ziek for new treatment
  P_nt[1, 1] <- 1 - P_nt[1, 2] - P_cp[1, 3]  # Adjust for remaining probability
  
  # Initialize arrays to store state proportions over time
  state_matrix_cp <- matrix(0, nrow = n_t + 1, ncol = n_states)
  state_matrix_nt <- matrix(0, nrow = n_t + 1, ncol = n_states)
  state_matrix_cp[1, ] <- c(1, 0, 0)  # Initial state: All individuals are 'Gezond'
  state_matrix_nt[1, ] <- c(1, 0, 0)
  
  # Propagate the state transitions over the time horizon
  for (t in 1:n_t) {
    state_matrix_cp[t + 1, ] <- state_matrix_cp[t, ] %*% P_cp
    state_matrix_nt[t + 1, ] <- state_matrix_nt[t, ] %*% P_nt
  }
  
  # Define utility and cost matrices
  utility_matrix <- matrix(c(u_gezond, u_ziek, u_dood), nrow = n_states, byrow = TRUE)
  cost_matrix <- matrix(c(c_gezond, c_ziek, c_dood), nrow = n_states, byrow = TRUE)
  
  # Calculate total costs and QALYs for both strategies using the utility and cost matrices
  total_QALYs_cp <- sum(state_matrix_cp %*% utility_matrix)
  total_QALYs_nt <- sum(state_matrix_nt %*% utility_matrix)
  total_cost_cp <- sum(state_matrix_cp %*% cost_matrix)
  total_cost_nt <- sum(state_matrix_nt %*% cost_matrix)
  
  # Return the total costs and QALYs for both strategies
  return(c(total_cost_cp, total_cost_nt, total_QALYs_cp, total_QALYs_nt))
}

