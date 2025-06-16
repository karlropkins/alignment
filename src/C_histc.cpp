#include <Rcpp.h>
using namespace Rcpp;

//   about Rcpp
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/

// histc <- function(values, edges) {
//   ledges <- length(edges)
//  nlow <-  length(values[values < edges[1]])
//  nhigh <- length(values[values > edges[ledges]])
//  bin <- NULL
//  for (i in seq_along(edges)[-ledges] ) {
//    binValues <- values[values >= edges[i] & values < edges[i+1]]
//    position <- match(binValues, values)
//    bin[position] <- i
//  }
//  bin <- bin[!is.na(bin)]
//  upper <- match(edges[ledges], values)
//    if(!is.na(upper)){
//      bin[length(bin) +1 ] = i+1
//    }
//    if( nlow > 0 ){
//      bin <- c(rep(0, nlow), bin)
//    }
//    if( nhigh > 0 ){
//      bin <- c(bin, rep(0, nhigh))
//    }
//    return(bin)
//}


// [[Rcpp::export]]
NumericVector C_histc(NumericVector values, NumericVector edges) {
  
  NumericVector nlow = values[values < edges[0]];
  NumericVector nhigh = values[values > edges[edges.size() - 1]];
  NumericVector bin(values.size(), NumericVector::get_na()); // get_na from Rccp 
  for (int i = 0; i < edges.size()-1; ++i) {
    for (int j = 0; j < values.size(); ++j) {
      if (values[j] >= edges[i] && values[j] < edges[i + 1]) {
        bin[j] = i + 1; // 1 + c++ bin index (ML>R>C inefficiency??)
      }
    }
  }
  bin = na_omit(bin);
  NumericVector upper(values.size(), NumericVector::get_na()); // from 
  for (int j = 0; j < values.size(); ++j) {
    if(values[j] == edges[edges.size() - 1]){
      upper[j] = j;
    }
  }  
  upper = na_omit(upper);
  if(upper.size() > 0) {
    bin.push_back(bin[bin.length()-1]+1);
  }
  if (nlow.size() > 0) {
    for (int k = 0; k < nlow.size(); ++k) {
      bin.push_front(0);
    }
  }
  if (nhigh.size() > 0) {
    for (int l = 0; l < nhigh.size(); ++l) {
      bin.push_back(0);
    }
  }
  return bin;
}

// C version of histc... in R  
// could improve/tidy this a lot but it works!!!!
// hurts my head thinking about this...


// this this should generate 0:9
// /*** R
// C_histc(1:10, 2:12)
// */
