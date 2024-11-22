#' @title State-Transition Model (without df conversion for mapply) 
#' @description Calculates state transitions, QALYs, and costs
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
#' f_model_b(df_input[1, ])
f_model_b <- function(params) {
  
  with(as.list(params), {
    
    # Initialize transition probability matrices
    a_P <- array(
      data = 0,
      dim = c(n_treatments, n_t, n_states, n_states),
      dimnames = list(v_treatments, 1:n_t, v_states, v_states)
    )
    
    # Transition probabilities for treatment 1
    # Transitions from "Gezond"
    a_P[1,, v_states[1], v_states[2]] <- p_gezond_ziek
    a_P[1,, v_states[1], v_states[3]] <- p_gezond_dood
    a_P[1,, v_states[1], v_states[1]] <- 1 - p_gezond_ziek - p_gezond_dood
    
    # Transitions from "Ziek"
    a_P[1,, v_states[2], v_states[1]] <- p_ziek_gezond
    a_P[1,, v_states[2], v_states[3]] <- p_ziek_dood
    a_P[1,, v_states[2], v_states[2]] <- 1 - p_ziek_gezond - p_ziek_dood
    
    # Transitions from "Dood"
    a_P[1,, v_states[3], v_states[3]] <- 1  # "Dood" is absorbing
    
    # Transition probabilities for treatment 2 (copy from treatment 1, adjust for relative risk)
    a_P[2,,,] <- a_P[1,,,]
    a_P[2,, v_states[1], v_states[2]] <- p_gezond_ziek * rr_gezond_ziek_t2_t1
    a_P[2,, v_states[1], v_states[3]] <- p_gezond_dood
    a_P[2,, v_states[1], v_states[1]] <- 1 - (p_gezond_ziek * rr_gezond_ziek_t2_t1) - p_gezond_dood
    
    # Initialize Markov trace
    a_TR <- array(
      data = NA,
      dim = c(n_treatments, n_t + 1, n_states),
      dimnames = list(v_treatments, 0:n_t, v_states)
    )
    a_TR[1, 1,] <- a_TR[2, 1,] <- c(1, 0, 0)  # Starting state: "Gezond"
    
    # State transitions using nested loops
    for (i_treatment in 1:n_treatments) {
      for (t in 1:n_t) {
        a_TR[i_treatment, t + 1, ] <- a_TR[i_treatment, t, ] %*% a_P[i_treatment, t, , ]
      }
    }
    
    # Create utility matrix
    m_u <- matrix(
      data = c(rep(u_gezond, n_t + 1),  # Utility for "Gezond"
               rep(u_ziek, n_t + 1),    # Utility for "Ziek"
               rep(u_dood, n_t + 1)),   # Utility for "Dood"
      nrow = n_t + 1,
      ncol = length(v_states)
    )
    
    # Create cost matrix
    m_c <- matrix(
      data = c(rep(c_gezond, n_t + 1),  # Costs for "Gezond"
               rep(c_ziek, n_t + 1),    # Costs for "Ziek"
               rep(c_dood, n_t + 1)),   # Costs for "Dood"
      nrow = n_t + 1,
      ncol = length(v_states)
    )
    
    # Estimate QALYs and costs
    v_E_t1 <- mapply(FUN = '%*%', t(a_TR[1, , ]), t(m_u))  # QALYs for treatment 1
    v_E_t2 <- mapply(FUN = '%*%', t(a_TR[2, , ]), t(m_u))  # QALYs for treatment 2
    v_C_t1 <- mapply(FUN = '%*%', t(a_TR[1, , ]), t(m_c))  # Costs for treatment 1
    v_C_t2 <- mapply(FUN = '%*%', t(a_TR[2, , ]), t(m_c))  # Costs for treatment 2
    
    # Return summed results
    return(c(sum(v_C_t1), sum(v_C_t2), sum(v_E_t1), sum(v_E_t2)))
    
  })
}
