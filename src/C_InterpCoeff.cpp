// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
using namespace Rcpp;
//   about Rcpp

// C version of InterCoeff... in R  
// could improve/tidy this a lot but *think* it works!!!!

//InterpCoeff <- function(n, nprime, offs, rtn) {
//  p <- length(nprime)
//  q <- n - 1
//  coeff <- matrix(nrow = p, ncol =  n)
//  index <- matrix(nrow = p, ncol =  n)
//  for (i in seq(1, p)) {
//    pp <- seq(1, nprime[i])
//    p <- seq(0, q) * (nprime[i] - 1)/q +1
//    k <- histc(p, pp)
//    k[p >= nprime[i] ] <- nprime[i] -1
//    coeff[i, ] <- ( p - pp[k] )
//    index[i, ] <- k -offs[i]
//  }
//  switch(rtn,
//         'coeff' = return(coeff),
//         'index' = return(index))
//}


// [[Rcpp::export]]
NumericMatrix c_InterpCoeff(NumericVector n, NumericVector nprime, NumericVector offs, CharacterVector rtn) {
  
  Environment pkg = Environment::namespace_env("alignment");
  Function c_histc = pkg["c_histc"];
  int p = nprime.length();
  int q = n[0] - 1;
  NumericMatrix coeff( p, n[0] ); 
  NumericMatrix index( p, n[0] );

  for(int i = 0; i < p; ++i) {
    IntegerVector pp = seq(1, nprime[i]);
    IntegerVector p2 = seq(0, q) ;
    NumericVector p3 = as<NumericVector>(p2) * (nprime[i] - 1)/q + 1; 
    NumericVector k  = c_histc(p3, pp);
    for(int j = 0; j < k.size(); ++j) {
      if(p3[j] >= nprime[i]){
        k[j] =nprime[i]-1;
      }
    }
    // Rcout << k << "    "; 
    for(int l = 0; l < k.size(); ++l) {
      coeff(i,l) = ( p3[l] - pp[k[l]-1] );
      index(i, l) = ( k[l] - offs[i] );
    }
  }
    
  if( rtn[0] == "coeff" ){
    return coeff;
  } else {
    return index;
  }

}

// /*** R
// this this should generate 
//      [,1] [,2] [,3]
// [1,]    0  0.5    1
// [2,]    0  0.0    1
// [3,]    0  0.5    1
// C_InterpCoeff(3, 2:4, -1:1, "coeff")
// */
