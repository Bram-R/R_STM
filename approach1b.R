# 1b Run model calculations with for loop 
for(x in 1:n_sim){ # loop to run model for each row of PSA inputs
  m_out_1b[x, ] <- f_model(x)
} # close for loop