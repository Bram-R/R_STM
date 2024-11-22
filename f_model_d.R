#' @title Optimized State-Transition Model (without df conversion, without mapply, without with(), etc) 
#' @description Calculates state transitions, QALYs, and costs
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
#' f_model_d(df_input[1, ])
f_model_d <- function(params) {
  
  # Initialize transition probability matrices
  a_transition <- array(
    data = 0,
    dim = c(n_treatments, n_t, n_states, n_states),
    dimnames = list(v_treatments, 1:n_t, v_states, v_states)
  )
  
  # Transition probabilities for treatment 1
  a_transition[1,, v_states[1], ] <- matrix(c(           # From health state 1 "Healthy"
    1 - params$p_healthy_sick - params$p_healthy_death,  # Stay in "Healthy"
    params$p_healthy_sick,                               # Transition to "Sick"
    params$p_healthy_death                               # Transition to "Death"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  a_transition[1,, v_states[2], ] <- matrix(c(           # From health state 2 "Sick"
    params$p_sick_healthy,                               # Transition to "Healthy"
    1 - params$p_sick_healthy - params$p_sick_death,     # Stay in "Sick"
    params$p_sick_death                                  # Transition to "Death"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  a_transition[1,, v_states[3], v_states[3]] <- 1        # "Death" is absorbing
  
  # Transition probabilities for treatment 2
  a_transition[2,,,] <- a_transition[1,,,]               # Copy from treatment 1
  a_transition[2,, v_states[1], ] <- matrix(c(                                         # From health state 1 "Healthy"
    1 - params$p_healthy_sick * params$rr_healthy_sick_t2_t1 - params$p_healthy_death, # Stay in "Healthy"
    params$p_healthy_sick * params$rr_healthy_sick_t2_t1,                              # Transition to "Sick"
    params$p_healthy_death                                                             # Transition to "Death"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  # Initialize Markov trace
  a_state_trace <- array(
    data = NA,
    dim = c(n_treatments, n_t + 1, n_states),
    dimnames = list(v_treatments, 0:n_t, v_states)
  )
  a_state_trace[1, 1,] <- a_state_trace[2, 1,] <- c(1, 0, 0)  # Starting state: "Healthy"
  
  # State transitions using a loop
  for (t in 1:n_t) {
    a_state_trace[1, t + 1, ] <- a_state_trace[1, t, ] %*% a_transition[1, t, , ]  # Treatment 1
    a_state_trace[2, t + 1, ] <- a_state_trace[2, t, ] %*% a_transition[2, t, , ]  # Treatment 2
  }
  
  # Utility and cost matrices
  m_utility <- matrix(c(params$u_healthy, # Utility for "Healthy"
                        params$u_sick,    # Utility for "Sick"
                        params$u_death),  # Utility for "Death"
                      nrow = n_t + 1, ncol = n_states, byrow = TRUE)
  m_cost <- matrix(c(params$c_healthy,    # Costs for "Healthy"
                     params$c_sick,       # Costs for "Sick"
                     params$c_death),     # Costs for "Death"
                   nrow = n_t + 1, ncol = n_states, byrow = TRUE)
  
  # Calculate QALYs and costs
  v_qalys <- rowSums(a_state_trace[, , ]  * rep(m_utility, each = n_treatments))
  v_costs <- rowSums(a_state_trace[, , ]  * rep(m_cost, each = n_treatments))
  
  # Return  results
  return(c(v_costs, v_qalys))
}
