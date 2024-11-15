#### Model calculations without df conversion for mapply ####

f_model_b <- function(params) {
  with(as.list(params), {
    
    ## set up transition probability matrix 
    # array with 4 dimensions (allowing transition matrix to depend on treatment and cycle)
    a_P <- array( 
      data = 0,
      dim = c(n_treatments, n_t, n_states, n_states), 
      dimnames = list(v_treatments, 1:n_t, v_states, v_states)  
    ) 
    
    # transitions from "Gezond"
    a_P[1,, v_states[1], v_states[2]] <- p_gezond_ziek
    a_P[1,, v_states[1], v_states[3]] <- p_gezond_dood
    a_P[1,, v_states[1], v_states[1]] <- 1 - p_gezond_ziek - p_gezond_dood
    
    # transitions from "Ziek"
    a_P[1,, v_states[2], v_states[1]] <- p_ziek_gezond
    a_P[1,, v_states[2], v_states[3]] <- p_ziek_dood
    a_P[1,, v_states[2], v_states[2]] <- 1 - p_ziek_gezond - p_ziek_dood
    
    # transitions from "Dood"
    a_P[1,, v_states[3], v_states[3]] <- 1 
    
    # copy transition matrix for new treatment
    a_P[2,,,] <- a_P[1,,,] 
    
    # transitions from "Gezond" for new treatment
    a_P[2,, v_states[1], v_states[2]] <- p_gezond_ziek * rr_gezond_ziek_t2_t1
    a_P[2,, v_states[1], v_states[3]] <- p_gezond_dood
    a_P[2,, v_states[1], v_states[1]] <- 1 - (p_gezond_ziek * rr_gezond_ziek_t2_t1) - p_gezond_dood
    
    ## state-transition model 
    a_TR <- array( 
      data = NA, 
      dim = c(n_treatments, n_t + 1, n_states), 
      dimnames = list(v_treatments, 0:n_t, v_states) 
    ) 
    
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
    
    m_u <- cbind( 
      rep(x = u_gezond, times = n_t + 1),
      rep(x = u_ziek, times = n_t + 1),
      rep(x = u_dood, times = n_t + 1)
    )
    
    # create cost matrix 
    m_c <- matrix( 
      data = 0, 
      nrow = n_t + 1, 
      ncol = length(v_states)
    ) # end matrix
    
    m_c <- cbind(
      rep(x = c_gezond, times = n_t + 1),
      rep(x = c_ziek, times = n_t + 1),
      rep(x = c_dood, times = n_t + 1)
    )
    
    # estimate QALYs and costs
    v_E_t1 <- mapply(FUN = '%*%', t(a_TR[1, , ]), t(m_u))
    v_E_t2 <- mapply(FUN = '%*%', t(a_TR[2, , ]), t(m_u))
    v_C_t1 <- mapply(FUN = '%*%', t(a_TR[1, , ]), t(m_c))
    v_C_t2 <- mapply(FUN = '%*%', t(a_TR[2, , ]), t(m_c))
    
    return(c(sum(v_C_t1) , sum(v_C_t2) , sum(v_E_t1), sum(v_E_t2))) # return model results
    
  }) # end with function  
} # end function
