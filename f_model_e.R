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
#' @return A numeric vector with the total costs and QALYs for each treatment: 
#'   \code{c(cost_t1, cost_t2, QALY_t1, QALY_t2)}.
#' @examples
#' f_model_e(df_input[1, ])
f_model_e <- function(params) {
  # Validate Rcpp source for parallel functionality
  # Rcpp::sourceCpp("f_propagate_states.cpp")
  
  # Initialize transition probability matrices
  a_P <- array(
    data = 0,
    dim = c(n_treatments, n_t, n_states, n_states),
    dimnames = list(v_treatments, 1:n_t, v_states, v_states)
  )
  
  # Transition probabilities for treatment 1
  a_P[1,, v_states[1], ] <- matrix(c(
    1 - params$p_gezond_ziek - params$p_gezond_dood,  # Stay in "Gezond"
    params$p_gezond_ziek,                             # Transition to "Ziek"
    params$p_gezond_dood                              # Transition to "Dood"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  a_P[1,, v_states[2], ] <- matrix(c(
    params$p_ziek_gezond,                             # Transition to "Gezond"
    1 - params$p_ziek_gezond - params$p_ziek_dood,    # Stay in "Ziek"
    params$p_ziek_dood                                # Transition to "Dood"
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  a_P[1,, v_states[3], v_states[3]] <- 1              # Remain in "Dood"
  
  # Transition probabilities for treatment 2
  a_P[2,,,] <- a_P[1,,,]
  a_P[2,, v_states[1], ] <- matrix(c(
    1 - params$p_gezond_ziek * params$rr_gezond_ziek_t2_t1 - params$p_gezond_dood,
    params$p_gezond_ziek * params$rr_gezond_ziek_t2_t1,
    params$p_gezond_dood
  ), nrow = n_t, ncol = n_states, byrow = TRUE)
  
  # Initialize Markov trace
  a_TR <- array(
    data = NA,
    dim = c(n_treatments, n_t + 1, n_states),
    dimnames = list(v_treatments, 0:n_t, v_states)
  )
  a_TR[1, 1,] <- a_TR[2, 1,] <- c(1, 0, 0)  # Starting state: "Gezond"
  
  # State transition using a Rcpp function
  a_TR[1, , ] <- f_propagate_states(a_TR[1, , ], a_P[1, , , ])
  a_TR[2, , ] <- f_propagate_states(a_TR[2, , ], a_P[2, , , ])

  # the f_propagate_states() calculations (in C++) are identical to
  # vectorized state transition using matrix multiplication
  # for (t in 1:n_t){ # loop through the number of cycles
  #   a_TR[1, t + 1, ] <- a_TR[1, t, ] %*% a_P[1, t, , ] # estimate next cycle (t + 1) of Markov trace for t1
  #   a_TR[2, t + 1, ] <- a_TR[2, t, ] %*% a_P[2, t, , ] # estimate next cycle (t + 1) of Markov trace for t2
  # } # close for loop for cycles
  
  # Utility and cost matrices
  m_u <- matrix(c(params$u_gezond, params$u_ziek, params$u_dood),
                nrow = n_t + 1, ncol = n_states, byrow = TRUE)
  m_c <- matrix(c(params$c_gezond, params$c_ziek, params$c_dood),
                nrow = n_t + 1, ncol = n_states, byrow = TRUE)
  
  # Calculate QALYs and costs
  m_res <- matrix(c(
    rowSums(a_TR[1, , ] * m_u),  # QALYs for treatment 1
    rowSums(a_TR[2, , ] * m_u),  # QALYs for treatment 2
    rowSums(a_TR[1, , ] * m_c),  # Costs for treatment 1
    rowSums(a_TR[2, , ] * m_c)   # Costs for treatment 2
  ), ncol = 4, byrow = FALSE)
  
  # Return aggregated results
  return(c(sum(m_res[, 3]), sum(m_res[, 4]), sum(m_res[, 1]), sum(m_res[, 2])))
}
