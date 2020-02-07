// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// Fragment Length Distribution Sampling

//' @export
// [[Rcpp::export]]
arma::mat fragLens_dist(int total_reads,
                        double mu, // Mean of normal dist
                        double SD, // SD of normal dist
                        int lenMin, // Minimum fragment length to consider
                        int lenMax){ // Maximum fragment length to consider
  
  int len = 0;
  int f = 0;
  int frag_round = 0;
  arma::mat dist_table((lenMax - lenMin + 1),2);
  arma::colvec lengths((lenMax - lenMin + 1));
  arma::colvec freq((lenMax - lenMin + 1));
  double frag = 0;
  arma::uvec index(1);
  
  for(len = 0; len < (lenMax - lenMin + 1); len++){
    lengths(len) = lenMin + len;
  }
  
  freq.zeros();
  
  for(f = 0; f < total_reads; f++){
    frag = R::rnorm(mu, SD);
    frag_round = round(frag);
    if(frag_round < lenMin){
      frag_round = lenMin;
    }else if(frag_round > lenMax){
      frag_round = lenMax;
    }
    
    index = find(lengths == frag_round); 
    freq(index) = freq(index) + 1.0;
    
  }
  
  dist_table.col(0) = freq;
  dist_table.col(1) = lengths;
  
  return(dist_table);
  
} // End fragLens_dist
  
// The End