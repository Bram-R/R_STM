#### Model calculations optimised (without df conversion, without mapply, without with(), etc) ####

f_model_d <- function(params) {
  
  ## set up transition probability matrix 
  # initialize array with 4 dimensions (allowing transition matrix to depend on treatment and cycle)
  a_P <- array( 
    data = 0,
    dim = c(n_treatments, n_t, n_states, n_states), 
    dimnames = list(v_treatments, 1:n_t, v_states, v_states)  
  ) 
  
  # vectorized assignment for transitions from "Gezond"
  a_P[1,, v_states[1],] <- matrix(c(1 - params$p_gezond_ziek - params$p_gezond_dood, # to "Gezond"
                                    params$p_gezond_ziek, # to "Ziek"
                                    params$p_gezond_dood), # to "Dood
                                  nrow = n_t, ncol = n_states, byrow = TRUE)
  
  # vectorized assignment for transitions from "Ziek"
  a_P[1,, v_states[2],] <- matrix(c(params$p_ziek_gezond, # to "Gezond"
                                    1 - params$p_ziek_gezond - params$p_ziek_dood, # to "Ziek"
                                    params$p_ziek_dood), # to "Dood
                                  nrow = n_t, ncol = n_states, byrow = TRUE)
  
  # transitions from "Dood"
  a_P[1,, v_states[3], v_states[3]] <- 1 
  
  # copy transition matrix for new treatment
  a_P[2,,,] <- a_P[1,,,] 
  
  # vectorized assignment for transitions from "Gezond" for new treatment
  a_P[2,, v_states[1],] <- matrix(c(1 - params$p_gezond_ziek * params$rr_gezond_ziek_t2_t1 - params$p_gezond_dood, # to "Gezond"
                                    params$p_gezond_ziek * params$rr_gezond_ziek_t2_t1, # to "Ziek"
                                    params$p_gezond_dood), # to "Dood
                                  nrow = n_t, ncol = n_states, byrow = TRUE)
  
  ## state-transition model 
  # initialize array with 3 dimensions 
  a_TR <- array( 
    data = NA, 
    dim = c(n_treatments, n_t + 1, n_states), 
    dimnames = list(v_treatments, 0:n_t, v_states) 
  ) 
  
  a_TR[1, 1,] <- a_TR[2, 1,] <- c(1, 0, 0) # set "Gezond" as starting health state
  
  # vectorized state transition using matrix multiplication
  for (t in 1:n_t){ # loop through the number of cycles
    a_TR[1, t + 1, ] <- a_TR[1, t, ] %*% a_P[1, t, , ] # estimate next cycle (t + 1) of Markov trace for t1
    a_TR[2, t + 1, ] <- a_TR[2, t, ] %*% a_P[2, t, , ] # estimate next cycle (t + 1) of Markov trace for t2
  } # close for loop for cycles
  
  ## calculate output
  # create utility matrix
  m_u <- matrix(c(params$u_gezond, 
                  params$u_ziek, 
                  params$u_dood), 
                nrow = n_t + 1, ncol = length(v_states), byrow = TRUE)
  
  # create cost matrix 
  m_c <- matrix(c(params$c_gezond, 
                  params$c_ziek, 
                  params$c_dood), 
                nrow = n_t + 1, ncol = length(v_states), byrow = TRUE)
  
  # estimate QALYs and costs with direct multiplication
  m_res <- matrix(c(rowSums(a_TR[1, , ] * m_u),  # QALYs for treatment 1
                    rowSums(a_TR[2, , ] * m_u),  # QALYs for treatment 2
                    rowSums(a_TR[1, , ] * m_c),  # Costs for treatment 1
                    rowSums(a_TR[2, , ] * m_c)), # Costs for treatment 2
                  ncol = 4, byrow = FALSE)
  
  return(c(sum(m_res[,3]) , sum(m_res[,4]) , sum(m_res[,1]), sum(m_res[,2]))) # return model results
  
} # end function
