f_model_a <- function(params) {
  #' @title State-Transition Model (with df conversion for mapply) 
  #' @description Calculates state transitions, QALYs, and costs
  #' @param params A named vector containing the following parameters:
  #'   \itemize{
  #'     \item \code{p_healthy_sick}: Transition probability from "Healthy" to "Sick".
  #'     \item \code{p_healthy_death}: Transition probability from "Healthy" to "Death".
  #'     \item \code{p_sick_healthy}: Transition probability from "Sick" to "Healthy".
  #'     \item \code{p_sick_death}: Transition probability from "Sick" to "Death".
  #'     \item \code{rr_healthy_sick_t2_t1}: Relative risk for new treatment.
  #'     \item \code{u_healthy}, \code{u_sick}, \code{u_death}: Utilities for each health state.
  #'     \item \code{c_healthy}, \code{c_sick}, \code{c_death}: Costs for each health state.
  #'   }
  #' @return A named vector with the total costs and QALYs for each treatment.
  #' @examples
  #' params <- setNames(c(0.33, 0.05, 0.20, 0.33, 0.82, 0.43, 0.38, 0.00, 7500, 10000, 0),
  #' c("p_healthy_sick", "p_healthy_death", "p_sick_healthy", "p_sick_death", "rr_healthy_sick_t2_t1",
  #'   "u_healthy", "u_sick", "u_death", "c_healthy", "c_sick", "c_death"))
  #' f_model_a(params)
    
  with(as.list(params), {
    
    # Initialize transition probability matrices
    a_transition <- array(
      data = 0,
      dim = c(n_treatments, n_t, n_states, n_states),
      dimnames = list(v_treatments, 1:n_t, v_states, v_states)
    )
    
    # Transition probabilities for treatment 1
    # From health state 1 "Healthy"
    a_transition[1,, v_states[1], v_states[1]] <- 1 - p_healthy_sick - p_healthy_death # Stay in health state 1 "Healthy"
    a_transition[1,, v_states[1], v_states[2]] <- p_healthy_sick                       # Transition to health state 2 "Sick"
    a_transition[1,, v_states[1], v_states[3]] <- p_healthy_death                      # Transition to health state 3 "Death"
    
    # From health state 2 "Sick"
    a_transition[1,, v_states[2], v_states[1]] <- p_sick_healthy                       # Transition to health state 1 "Healthy"
    a_transition[1,, v_states[2], v_states[2]] <- 1 - p_sick_healthy - p_sick_death    # Stay in health state 2 "Sick"
    a_transition[1,, v_states[2], v_states[3]] <- p_sick_death                         # Transition to health state 3 "Death"
    
    
    # From health state 3 "Death"
    a_transition[1,, v_states[3], v_states[3]] <- 1                                    # Heath state 3 "Death" is absorbing
    
    # Transition probabilities for treatment 2
    # Copy from treatment 1
    a_transition[2,,,] <- a_transition[1,,,]
    
    # From health state 1 "Healthy"
    a_transition[2,, v_states[1], v_states[1]] <- 1 - (p_healthy_sick * rr_healthy_sick_t2_t1) - p_healthy_death  # Stay in health state 1 "Healthy"
    a_transition[2,, v_states[1], v_states[2]] <- p_healthy_sick * rr_healthy_sick_t2_t1                          # Transition to health state 2 "Sick"
    a_transition[2,, v_states[1], v_states[3]] <- p_healthy_death                                                 # Transition to health state 3 "Death"
    
    # Initialize Markov trace
    a_state_trace <- array(
      data = NA,
      dim = c(n_treatments, n_t + 1, n_states),
      dimnames = list(v_treatments, 0:n_t, v_states)
    )
    a_state_trace[1, 1,] <- a_state_trace[2, 1,] <- c(1, 0, 0) # Starting health state: 1 "Healthy"
    
    # State transitions using nested loops
    for (i_treatment in 1:n_treatments) {
      for (t in 1:n_t) {
        a_state_trace[i_treatment, t + 1, ] <- a_state_trace[i_treatment, t, ] %*% a_transition[i_treatment, t, , ]
      }
    }
    
    # Create utility matrix
    m_utility <- matrix( 
      data = NA, 
      nrow = n_t + 1, 
      ncol = length(v_states)
    ) 
    m_utility <- cbind( 
      rep(x = u_healthy, times = n_t + 1), # Utility for health state 1 "Healthy"
      rep(x = u_sick, times = n_t + 1),    # Utility for health state 2 "Sick"
      rep(x = u_death, times = n_t + 1)    # Utility for health state 3 "Death"
    )
    
    # Create cost matrix 
    m_cost <- matrix( 
      data = 0, 
      nrow = n_t + 1, 
      ncol = length(v_states)
    ) 
    m_cost <- cbind(
      rep(x = c_healthy, times = n_t + 1), # Costs for health state 1 "Healthy"
      rep(x = c_sick, times = n_t + 1),    # Costs for health state 2 "Sick"
      rep(x = c_death, times = n_t + 1)    # Costs for health state 3 "Death"
    )
    
    # Calculate QALYs and costs
    v_E_t1 <- mapply(FUN = '%*%', as.data.frame(t(a_state_trace[1, , ])), as.data.frame(t(m_utility))) # QALYs for treatment 1
    v_E_t2 <- mapply(FUN = '%*%', as.data.frame(t(a_state_trace[2, , ])), as.data.frame(t(m_utility))) # QALYs for treatment 2
    v_C_t1 <- mapply(FUN = '%*%', as.data.frame(t(a_state_trace[1, , ])), as.data.frame(t(m_cost))) # Costs for treatment 1
    v_C_t2 <- mapply(FUN = '%*%', as.data.frame(t(a_state_trace[2, , ])), as.data.frame(t(m_cost))) # Costs for treatment 2
    
    # Named vector with results
    v_results <- setNames(c(sum(v_C_t1), sum(v_C_t2), sum(v_E_t1), sum(v_E_t2)),            # Combined results
                          c(paste0("Cost_", v_treatments), paste0("QALY_", v_treatments)))  # Result names 
    
    return(v_results)
  })
}

# Test
# # Defined names vector for input parameters
# v_params_test <- c(p_healthy_sick = 0.4,
#                    p_healthy_death = 0.1,
#                    p_sick_healthy = 0.2,
#                    p_sick_death = 0.3,
#                    rr_healthy_sick_t2_t1 = 0.5,
#                    u_healthy = 1,
#                    u_sick = 0.5,
#                    u_death = 0,
#                    c_healthy = 100,
#                    c_sick = 1000,
#                    c_death = 0)
# 
# # Test the model
# f_model_a(v_params_test)
