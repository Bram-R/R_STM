#' @title Optimized State-Transition Model using Rcpp
#' @description Calculates state transitions, QALYs, and costs using Rcpp for efficiency.
#' @param params A data frame containing the following parameters:
#'   \itemize{
#'     \item \code{p_healthy_sick}: Transition probability from "Healthy" to "Sick".
#'     \item \code{p_healthy_death}: Transition probability from "Healthy" to "Death".
#'     \item \code{p_sick_healthy}: Transition probability from "Sick" to "Healthy".
#'     \item \code{p_sick_death}: Transition probability from "Sick" to "Death".
#'     \item \code{rr_healthy_sick_t2_t1}: Relative risk for new treatment.
#'     \item \code{u_healthy}, \code{u_sick}, \code{u_death}: Utilities for each health state.
#'     \item \code{c_healthy}, \code{c_sick}, \code{c_death}: Costs for each health state.
#'   }
#' @return A numeric vector with the total costs and QALYs for each treatment.
#' @examples
#' f_model_e(df_input[1, ])
f_model_e <- function(params) {
  # Validate Rcpp source for parallel functionality
  # Rcpp::sourceCpp("f_propagate_states.cpp")
  
  # Initialize transition probability matrices
  a_transition <- array(
    data = 0,
    dim = c(n_treatments, n_t, n_states, n_states),
    dimnames = list(v_treatments, 1:n_t, v_states, v_states)
  )
  
  # Transition probabilities for treatment 1
  a_transition[1, , v_states[1], ] <- matrix(c(          # From health state 1 "Healthy"
    1 - params$p_healthy_sick - params$p_healthy_death,  # Stay in health state 1 "Healthy"
    params$p_healthy_sick,                               # Transition to health state 2 "Sick"
    params$p_healthy_death                               # Transition to health state 3 "Death"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  a_transition[1, , v_states[2], ] <- matrix(c(          # From health state 2 "Sick"
    params$p_sick_healthy,                               # Transition to health state 1"Healthy"
    1 - params$p_sick_healthy - params$p_sick_death,     # Stay in health state 2 "Sick"
    params$p_sick_death                                  # Transition to health state 3 "Death"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  a_transition[1, , v_states[3], v_states[3]] <- 1       # Heath state 3 "Death" is absorbing
  
  # Transition probabilities for treatment 2
  a_transition[2, , , ] <- a_transition[1, , , ]         # Copy from treatment 1
  a_transition[2, , v_states[1], ] <- matrix(c(                                        # From health state 1 "Healthy"
    1 - params$p_healthy_sick * params$rr_healthy_sick_t2_t1 - params$p_healthy_death, # Stay in health state 1 "Healthy"
    params$p_healthy_sick * params$rr_healthy_sick_t2_t1,                              # Transition to health state 2 "Sick"
    params$p_healthy_death                                                             # Transition to health state 3 "Death"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  # Initialize Markov trace
  a_state_trace <- array(
    data = NA,
    dim = c(n_treatments, n_t + 1, n_states),
    dimnames = list(v_treatments, 0:n_t, v_states)
  )
  a_state_trace[, 1, ] <- matrix(c(1, 0, 0), nrow = n_treatments, ncol = n_states, byrow = TRUE) # Starting health state: 1 "Healthy"
  
  # State transition using a Rcpp function with a loop for treatment
  for (i_treatment in 1:n_treatments) {
    a_state_trace[i_treatment, , ] <- f_propagate_states(a_state_trace[i_treatment, , ], a_transition[i_treatment, , , ])
  }
  
  # Utility and cost matrices
  m_utility <- matrix(c(params$u_healthy, # Utility for health state 1 "Healthy"
                        params$u_sick,    # Utility for health state 2 "Sick"
                        params$u_death),  # Utility for health state 3 "Death"
                      nrow = n_t + 1, ncol = n_states, byrow = TRUE)
  m_cost <- matrix(c(params$c_healthy,    # Costs for health state 1"Healthy"
                     params$c_sick,       # Costs for health state 2 "Sick"
                     params$c_death),     # Costs for health state 3 "Death"
                   nrow = n_t + 1, ncol = n_states, byrow = TRUE)
  
  # Calculate QALYs and costs
  v_qalys <- rowSums(a_state_trace[, , ]  * rep(m_utility, each = n_treatments))
  v_costs <- rowSums(a_state_trace[, , ]  * rep(m_cost, each = n_treatments))
  
  # Return  results
  return(c(v_costs, v_qalys))
}
