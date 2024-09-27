#### Model calculations without df conversion for mapply ####

f_model_matrix <- function(x) {

    ## set up transition probability matrix 
    # array with 4 dimensions (allowing transition matrix to depend on treatment and cycle)
    a_P <- array( # array with 4 dimensions (treatment, cycles, v_states x v_states tp matrix)
      data = 0,
      dim = c(n_treatments, n_t, n_states, n_states), 
      dimnames = list(v_treatments, 1:n_t, v_states, v_states)  
    ) # close array
    
    # transitions from "Gezond"
    a_P[1,, v_states[1], v_states[2]] <- x[1]
    a_P[1,, v_states[1], v_states[3]] <- x[2]
    a_P[1,, v_states[1], v_states[1]] <- 1 - x[1] - x[2]
    
    # transitions from "Ziek"
    a_P[1,, v_states[2], v_states[1]] <- x[3]
    a_P[1,, v_states[2], v_states[3]] <- x[4]
    a_P[1,, v_states[2], v_states[2]] <- 1 - x[3] - x[4]
    
    # transitions from "Dood"
    a_P[1,, v_states[3], v_states[3]] <- 1 
    
    # copy transition matrix for new treatment
    a_P[2,,,] <- a_P[1,,,] 
    
    # transitions from "Gezond" for new treatment
    a_P[2,, v_states[1], v_states[2]] <- x[1] * x[5]
    a_P[2,, v_states[1], v_states[3]] <- x[2]
    a_P[2,, v_states[1], v_states[1]] <- 1 - (x[1] * x[5]) - x[2]
    
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
      rep(x = x[6], times = n_t + 1),
      rep(x = x[7], times = n_t + 1),
      rep(x = x[8], times = n_t + 1)
    )
    
    # create cost matrix 
    m_c <- matrix( 
      data = 0, 
      nrow = n_t + 1, 
      ncol = length(v_states)
    ) # end matrix
    
    m_c <- cbind(
      rep(x = x[9], times = n_t + 1),
      rep(x = x[10], times = n_t + 1),
      rep(x = x[11], times = n_t + 1)
    )
    
    # estimate QALYs and costs
    v_E_t1 <- mapply(FUN = '%*%', t(a_TR[1, , ]), t(m_u))
    v_E_t2 <- mapply(FUN = '%*%', t(a_TR[2, , ]), t(m_u))
    v_C_t1 <- mapply(FUN = '%*%', t(a_TR[1, , ]), t(m_c))
    v_C_t2 <- mapply(FUN = '%*%', t(a_TR[2, , ]), t(m_c))
    
    return(c(sum(v_C_t1) , sum(v_C_t2) , sum(v_E_t1), sum(v_E_t2))) # return model results
    
} # end function
