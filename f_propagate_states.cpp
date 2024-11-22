#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Propagate States in a Markov Model
 //' @description This function propagates the state probabilities through a Markov model 
 //'   using a transition probability cube.
 //' @param a_TR An \code{arma::mat} containing the initial Markov trace with dimensions 
 //'   \code{(n_t + 1) x n_states}, where \code{n_t} is the number of time steps and 
 //'   \code{n_states} is the number of states.
 //' @param a_P An \code{arma::cube} containing transition probability matrices with 
 //'   dimensions \code{n_t x n_states x n_states}.
 //' @return An \code{arma::mat} containing the updated Markov trace after propagating 
 //'   through all time steps.
 //' @examples
 //' \dontrun{
 //' arma::mat trace = arma::zeros<arma::mat>(n_t + 1, n_states);
 //' trace.row(0) = initial_state; // Set initial state
 //' arma::cube transition_cube = ...; // Define transition probabilities
 //' arma::mat result = f_propagate_states(trace, transition_cube);
 //' }
 // [[Rcpp::export]]
 arma::mat f_propagate_states(arma::mat a_TR, const arma::cube& a_P) {
   
   // Number of time steps and states
   int n_t = a_P.n_rows;       // Number of time steps
   int n_states = a_TR.n_cols; // Number of states
   
   // Loop through each time step
   for (int t = 0; t < n_t; ++t) {
     
     // Extract the transition matrix for the current time step
     arma::mat transition_matrix = a_P.tube(t, 0, t, n_states - 1);
     
     // Update the state probabilities for the next time step
     a_TR.row(t + 1) = a_TR.row(t) * transition_matrix;
   }
   
   // Return the propagated Markov trace
   return a_TR;
 }
