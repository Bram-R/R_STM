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
#' f_model_c(df_input[1, ])

f_model_c <- function(params) {
  
  ## set up transition probability matrix 
  # array with 4 dimensions (allowing transition matrix to depend on treatment and cycle)
  a_P <- array( # array with 4 dimensions (treatment, cycles, v_states x v_states tp matrix)
    data = 0,
    dim = c(n_treatments, n_t, n_states, n_states), 
    dimnames = list(v_treatments, 1:n_t, v_states, v_states)  
  ) # close array
  
  # transitions from "Gezond"
  a_P[1,, v_states[1], v_states[2]] <- params[1]
  a_P[1,, v_states[1], v_states[3]] <- params[2]
  a_P[1,, v_states[1], v_states[1]] <- 1 - params[1] - params[2]
  
  # transitions from "Ziek"
  a_P[1,, v_states[2], v_states[1]] <- params[3]
  a_P[1,, v_states[2], v_states[3]] <- params[4]
  a_P[1,, v_states[2], v_states[2]] <- 1 - params[3] - params[4]
  
  # transitions from "Dood"
  a_P[1,, v_states[3], v_states[3]] <- 1 
  
  # copy transition matrix for new treatment
  a_P[2,,,] <- a_P[1,,,] 
  
  # transitions from "Gezond" for new treatment
  a_P[2,, v_states[1], v_states[2]] <- params[1] * params[5]
  a_P[2,, v_states[1], v_states[3]] <- params[2]
  a_P[2,, v_states[1], v_states[1]] <- 1 - (params[1] * params[5]) - params[2]
  
  ## state-transition model 
  a_TR <- array( # array with 3 dimensions (treatment, cycles, health states)
    data = NA, 
    dim = c(n_treatments, n_t + 1, n_states), 
    dimnames = list(v_treatments, 0:n_t, v_states) 
  ) # close array  
  
  a_TR[1, 1,] <- a_TR[2, 1,] <- c(1, 0, 0) # set "Gezond" as starting health state
  
  for (i_treatment in 1:n_treatments){ # loop through the treatment options
    for (t in 1:n_t){ # loop through the number of cycles
      a_TR[i_treatment, t + 1, ] <- a_TR[i_treatment, t, ] %*% a_P[i_treatment, t, , ] # estimate next cycle (t + 1) of Markov trace 
    } # close for loop for cycles
  } # close for loop for treatments
  
  ## calculate output
  # create utility matrix
  m_u <- matrix( 
    data = NA, 
    nrow = n_t + 1, 
    ncol = length(v_states)
  ) # end matrix
  
  # create utility matrix 
  m_u <- cbind( 
    rep(x = params[6], times = n_t + 1),
    rep(x = params[7], times = n_t + 1),
    rep(x = params[8], times = n_t + 1)
  )
  
  # create cost matrix 
  m_c <- matrix( 
    data = 0, 
    nrow = n_t + 1, 
    ncol = length(v_states)
  ) # end matrix
  
  m_c <- cbind(
    rep(x = params[9], times = n_t + 1),
    rep(x = params[10], times = n_t + 1),
    rep(x = params[11], times = n_t + 1)
  )
  
  # estimate QALYs and costs
  v_E_t1 <- mapply(FUN = '%*%', t(a_TR[1, , ]), t(m_u))
  v_E_t2 <- mapply(FUN = '%*%', t(a_TR[2, , ]), t(m_u))
  v_C_t1 <- mapply(FUN = '%*%', t(a_TR[1, , ]), t(m_c))
  v_C_t2 <- mapply(FUN = '%*%', t(a_TR[2, , ]), t(m_c))
  
  return(c(sum(v_C_t1) , sum(v_C_t2) , sum(v_E_t1), sum(v_E_t2))) # return model results
  
} # end function
