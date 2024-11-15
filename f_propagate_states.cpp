#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat f_propagate_states(arma::mat a_TR, const arma::cube& a_P) {
  int n_t = a_P.n_rows;    // Number of time steps (100)
  int n_states = a_TR.n_cols; // Number of states (3)
  
  // Iterate through each time step
  for (int t = 0; t < n_t; ++t) {
    // Extract the transition matrix for time step t (size 3x3)
    arma::mat transition_matrix = a_P.tube(t, 0, t, n_states - 1);
    // Perform matrix multiplication for the current time step
    a_TR.row(t + 1) = a_TR.row(t) * transition_matrix;
  }
  
  return a_TR;
}
