#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat f_propagate_states(const arma::mat& P, int n_t) {
  int n_states = P.n_rows;
  arma::mat state_matrix(n_t + 1, n_states, arma::fill::zeros);
  state_matrix(0, 0) = 1; // Initial state: all individuals in 'Gezond'
  
  // Propagate states over time
  for (int t = 0; t < n_t; ++t) {
    state_matrix.row(t + 1) = state_matrix.row(t) * P;
  }
  
  return state_matrix;
}
