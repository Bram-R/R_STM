#' @title State-Transition Model (without df conversion for mapply, without with()) 
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
#' @return A numeric vector with the total costs and QALYs for each treatment:
#'   \code{c(cost_t1, cost_t2, QALY_t1, QALY_t2)}.
#' @examples
#' f_model_c(as.matrix(df_input[1, ]))
f_model_c <- function(params) {
  
  # Initialize transition probability matrices
  a_transition <- array(
    data = 0,
    dim = c(n_treatments, n_t, n_states, n_states),
    dimnames = list(v_treatments, 1:n_t, v_states, v_states)
  )
  
  # Transition probabilities for treatment 1
  # Transitions from "Healthy"
  a_transition[1,, v_states[1], v_states[2]] <- params[1]  # Transition to "Sick"
  a_transition[1,, v_states[1], v_states[3]] <- params[2]  # Transition to "Death"
  a_transition[1,, v_states[1], v_states[1]] <- 1 - params[1] - params[2]  # Stay in "Healthy"
  
  # Transitions from "Sick"
  a_transition[1,, v_states[2], v_states[1]] <- params[3]  # Transition to "Healthy"
  a_transition[1,, v_states[2], v_states[3]] <- params[4]  # Transition to "Death"
  a_transition[1,, v_states[2], v_states[2]] <- 1 - params[3] - params[4]  # Stay in "Sick"
  
  # Transitions from "Death"
  a_transition[1,, v_states[3], v_states[3]] <- 1  # "Death" is absorbing
  
  # Transition probabilities for treatment 2
  a_transition[2,,,] <- a_transition[1,,,]
  a_transition[2,, v_states[1], v_states[2]] <- params[1] * params[5]  # Adjusted transition to "Sick"
  a_transition[2,, v_states[1], v_states[3]] <- params[2]  # Transition to "Death"
  a_transition[2,, v_states[1], v_states[1]] <- 1 - (params[1] * params[5]) - params[2]  # Stay in "Healthy"
  
  # Initialize Markov trace
  a_state_trace <- array(
    data = NA,
    dim = c(n_treatments, n_t + 1, n_states),
    dimnames = list(v_treatments, 0:n_t, v_states)
  )
  a_state_trace[1, 1,] <- a_state_trace[2, 1,] <- c(1, 0, 0)  # Starting state: "Healthy"
  
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
    rep(x = params[6], times = n_t + 1), # Utility for "Healthy"
    rep(x = params[7], times = n_t + 1), # Utility for "Sick"
    rep(x = params[8], times = n_t + 1)  # Utility for "Death"
  )
  
  # Create cost matrix 
  m_cost <- matrix( 
    data = 0, 
    nrow = n_t + 1, 
    ncol = length(v_states)
  ) 
  
  m_cost <- cbind(
    rep(x = params[9], times = n_t + 1),  # Costs for "Healthy"
    rep(x = params[10], times = n_t + 1), # Costs for "Sick"
    rep(x = params[11], times = n_t + 1)  # Costs for "Death"
  )
  
  # Estimate QALYs and costs
  v_E_t1 <- mapply(FUN = '%*%', t(a_state_trace[1, , ]), t(m_utility))  # QALYs for treatment 1
  v_E_t2 <- mapply(FUN = '%*%', t(a_state_trace[2, , ]), t(m_utility))  # QALYs for treatment 2
  v_C_t1 <- mapply(FUN = '%*%', t(a_state_trace[1, , ]), t(m_cost))     # Costs for treatment 1
  v_C_t2 <- mapply(FUN = '%*%', t(a_state_trace[2, , ]), t(m_cost))     # Costs for treatment 2
  
  # Return summed results
  return(c(sum(v_C_t1), sum(v_C_t2), sum(v_E_t1), sum(v_E_t2)))
}
