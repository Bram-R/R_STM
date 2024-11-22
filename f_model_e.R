#' @title Optimized State-Transition Model using Rcpp
#' @description Calculates state transitions, QALYs, and costs using Rcpp for efficiency.
#' @param params A data frame containing the following parameters:
#'   \itemize{
#'     \item \code{p_gezond_ziek}: Transition probability from "Gezond" to "Ziek".
#'     \item \code{p_gezond_dood}: Transition probability from "Gezond" to "Dood".
#'     \item \code{p_ziek_gezond}: Transition probability from "Ziek" to "Gezond".
#'     \item \code{p_ziek_dood}: Transition probability from "Ziek" to "Dood".
#'     \item \code{rr_gezond_ziek_t2_t1}: Relative risk for new treatment.
#'     \item \code{u_gezond}, \code{u_ziek}, \code{u_dood}: Utilities for each health state.
#'     \item \code{c_gezond}, \code{c_ziek}, \code{c_dood}: Costs for each health state.
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
  a_transition[1,, v_states[1], ] <- matrix(c(
    1 - params$p_gezond_ziek - params$p_gezond_dood,  # Stay in "Gezond"
    params$p_gezond_ziek,                             # Transition to "Ziek"
    params$p_gezond_dood                              # Transition to "Dood"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  a_transition[1,, v_states[2], ] <- matrix(c(
    params$p_ziek_gezond,                             # Transition to "Gezond"
    1 - params$p_ziek_gezond - params$p_ziek_dood,    # Stay in "Ziek"
    params$p_ziek_dood                                # Transition to "Dood"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  a_transition[1,, v_states[3], v_states[3]] <- 1  # "Dood" is absorbing
  
  # Transition probabilities for treatment 2
  a_transition[2,,,] <- a_transition[1,,,]
  a_transition[2,, v_states[1], ] <- matrix(c(
    1 - params$p_gezond_ziek * params$rr_gezond_ziek_t2_t1 - params$p_gezond_dood, # Stay in "Gezond"
    params$p_gezond_ziek * params$rr_gezond_ziek_t2_t1,                            # Transition to "Ziek"
    params$p_gezond_dood                                                           # Transition to "Dood"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  # Initialize Markov trace
  a_state_trace <- array(
    data = NA,
    dim = c(n_treatments, n_t + 1, n_states),
    dimnames = list(v_treatments, 0:n_t, v_states)
  )
  a_state_trace[1, 1,] <- a_state_trace[2, 1,] <- c(1, 0, 0)  # Starting state: "Gezond"
  
  # State transition using a Rcpp function
  a_state_trace[1, , ] <- f_propagate_states(a_state_trace[1, , ], a_transition[1, , , ])
  a_state_trace[2, , ] <- f_propagate_states(a_state_trace[2, , ], a_transition[2, , , ])
  
  # the f_propagate_states() calculations (in C++) are identical to
  # State transitions using a loop
  # for (t in 1:n_t) {
  #   a_state_trace[1, t + 1, ] <- a_state_trace[1, t, ] %*% a_transition[1, t, , ]  # Treatment 1
  #   a_state_trace[2, t + 1, ] <- a_state_trace[2, t, ] %*% a_transition[2, t, , ]  # Treatment 2
  # }
  
  # Utility and cost matrices
  m_utility <- matrix(c(params$u_gezond, # Utility for "Gezond"
                        params$u_ziek,   # Utility for "Ziek"
                        params$u_dood),  # Utility for "Dood"
                      nrow = n_t + 1, ncol = n_states, byrow = TRUE)
  m_cost <- matrix(c(params$c_gezond, # Costs for "Gezond"
                     params$c_ziek,   # Costs for "Ziek"
                     params$c_dood),  # Costs for "Dood"
                   nrow = n_t + 1, ncol = n_states, byrow = TRUE)
  
  # Calculate QALYs and costs
  v_qalys <- rowSums(a_state_trace[, , ]  * rep(m_utility, each = n_treatments))
  v_costs <- rowSums(a_state_trace[, , ]  * rep(m_cost, each = n_treatments))
  
  # Return  results
  return(c(v_costs, v_qalys))
}
