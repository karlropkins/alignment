#include <Rcpp.h>
using namespace Rcpp;

//   about Rcpp
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/

//cow as R////

// cow <- function(Ta, X, Seg, Slack, Options) {

// output list 
// https://stackoverflow.com/questions/72176779/how-to-return-objects-of-different-types-from-rcpparmadillo-function
// Rcpp::NumericMatrix  
// https://teuder.github.io/rcpp4everyone_en/100_matrix.html
// running R code from C
// https://teuder.github.io/rcpp4everyone_en/230_R_function.html 



// [[Rcpp::export]]
List C_cow(NumericVector Ta, NumericMatrix X, NumericVector Seg, NumericVector Slack, NumericVector Options) {
  
  // setup
  NumericVector dimX = X.attr("dim");
  Seg = floor(Seg);
  int npT = Ta.length();    // number of points in target
  NumericVector nSeg(1); 
  if(Seg.length()==1){
    // calculate number of segments  
    if(Options[2]==1){
      nSeg[0] = floor((npT - 1)/Seg[0]);
    } else {
      nSeg[0] = floor((npT - 1) / (Seg[0] - 1));
    }
  } else{
    nSeg[0] = Seg.length()-1;
  }
  NumericMatrix len_segs(2, nSeg[0]);
  // populate len_segs
  if (Seg.length()==1){ // Seg.length 1
    if(Options[2]==1){
      for(int i=0; i < nSeg[0]; ++i){
        len_segs(0,i) = ( floor((npT - 1)/ nSeg[0]) );
        len_segs(1,i) = ( floor((dimX[1] - 1)/ nSeg[0]) );
      }
      // from R*  print('Segment lengh adjusted to the best cover the remainders')
      Rcout << "Segment length adjusted to the best cover the remainders" << "\n"; 
    } else {
      for(int i=0; i < nSeg[0]; ++i){
        len_segs(0,i) = ( Seg[0] - 1 );
        len_segs(1,i) = ( Seg[0] - 1 );
      }
      // from R* if (floor((dimX[2] -1) / (Seg - 1)) != nSeg){
      if (floor((dimX[1]-1)/(Seg[0]-1)) != nSeg[0]){
        stop("For non-fixed segment lengths the target and the signal do not have the same number of segments. Try option 3 set to T");
      }
    }
    NumericVector tmp(1);
    tmp[0] = ( (npT-1) - ( floor((npT-1)/len_segs(0, 0)) *  len_segs(0, 0) ) ); 
    if( tmp[0] > 0) {
      len_segs(0, nSeg[0]-1) = len_segs(0, nSeg[0]-1) + tmp[0];
      if (Options[0] == 1) {
        Rcout << "Segments: " << len_segs(0,0) + 1 << " Points: " << nSeg[0]-1 << "\n";
      }
    } else {
      if (Options[0] == 1) {
        Rcout << "Segments: " << len_segs(1,0) + 1 << " Points: " << nSeg[0] << "\n";
      }
    }
    tmp[0] = ( (dimX[1]-1) - ( floor((dimX[1]-1)/len_segs(1, 0)) *  len_segs(1, 0) ) );
    if(tmp[0] > 0) {
      len_segs(1, nSeg[0]-1) = len_segs(1, nSeg[0]-1) + tmp[0];
    }
    
  } else {  // seg.length not 1
    if ( (Seg[0] != 1) | (Seg[Seg.length()-1] != npT) ) {
      stop("End points must be equal to 1 and to the length of the target"); 
    }
    for(int i=0; i < nSeg[0]; ++i){
      len_segs(0,i) = ( Seg[i+1] - Seg[i] );
      len_segs(1,i) = ( Seg[i+1] - Seg[i] );
    }
    // from R* if (!all(len_segs > 2)) {
    if(min(len_segs) < 2){
      stop("All segments must contain at least two points");
    }
  }
  
  NumericMatrix Xwarped(dimX[0], dimX[1]);  // Initialise matrix of warped signals
  NumericVector X0 = X(0,_);  // this is the first column of the supplied matrix
  
  NumericVector bT(nSeg[0]+1);
  bT[0] = 1;
  for(int i=1; i<nSeg[0]+1; ++i){
    bT[i] = bT[i-1] + len_segs(0,i-1);
  }
  NumericVector bP(nSeg[0]+1);
  bP[0] = 1;
  for(int i=1; i<nSeg[0]+1; ++i){
    bP[i] = bP[i-1] + len_segs(1,i-1);
  }
  NumericMatrix Warping( dimX[0], nSeg[0]+1) ; 
  std::fill(Warping.begin(), Warping.end(), NumericVector::get_na());
  
  // check slack
  if(Slack.length() > 1 ){
    if(Slack.length() <- nSeg[0] ){
      stop("The number of slack parameters is not equal to the number of optimised segments");
    }
    stop("Multiple slacks have not been implemented yet"); 
  }
  IntegerVector Slacks_vec = seq(-Slack[0], Slack[0]);
  
  
  //constraints
  IntegerVector offs_tmp = seq(0, nSeg[0]) * Slack[0]; 
  NumericMatrix offs(2, offs_tmp.size()); 
  for(int i=0; i<offs_tmp.length(); ++i){
    offs(0,i) = -offs_tmp[i]; 
    offs(1,i) = offs_tmp[i];
  }
  //I simplified next bit...
  //IntegerVector Bounds_ind = seq(0, nSeg[0] ); // don't need to use
  // TO LOOK AT 
  // Can I drop Bounds_a and Bounds_b and work striaght into Bounds
  //    would drop two matrices...
  NumericMatrix Bounds_a(2, offs_tmp.size()); 
  for(int i=0; i<offs_tmp.length(); ++i){
    Bounds_a(0,i) = bP[i] + offs(0,i);
    Bounds_a(1,i) = bP[i] + offs(1,i);      
  }
  NumericMatrix Bounds_b(2, offs_tmp.size()); 
  for(int i=0; i < offs_tmp.length(); ++i){
    Bounds_b(0,i) = bP[i] + offs(0,offs_tmp.size() -i -1);
    Bounds_b(1,i) = bP[i] + offs(1,offs_tmp.size() -i -1);      
  }
  NumericMatrix Bounds(2, offs_tmp.size()); 
  for(int i=0; i < offs_tmp.length(); ++i){
    if(Bounds_b(0,i) > Bounds_a(0,i)){
      Bounds(0,i) = Bounds_b(0,i);
    } else{
      Bounds(0,i) = Bounds_a(0,i);
    }
    if(Bounds_b(1,i) < Bounds_a(1,i)){
      Bounds(1,i) = Bounds_b(1,i);
    } else {
      Bounds(1,i) = Bounds_a(1,i);
    }
  }
  
  // ignoring option 4; coping daniels methods...
  //    could come back to this
  
  // calculate diffs
  NumericMatrix Xdiff(dimX[0], dimX[1]-1 ); 
  for(int i=0; i<dimX[0]; ++i){
    for(int j=1; j<dimX[1]; ++j){
      Xdiff(i,j-1) = X(i,j) - X(i,j-1);
    }
  }
  //  Calculate coefficients and indexes for interpolation ####
  List Int_Coeff(nSeg[0]); 
  List Int_Index(nSeg[0]);
  
  Environment pkg = Environment::namespace_env("alignment");
  Function C_InterpCoeff = pkg["C_InterpCoeff"];
  
  if(Seg.length()  == 1){
    // supplied seg length 1
    for(int i=0; i < nSeg[0]-1; ++i){
      Int_Coeff[i] = C_InterpCoeff(len_segs(0, 0) +1, len_segs(1, 0) + Slacks_vec + 1, Slacks_vec, "coeff");
      Int_Index[i] = C_InterpCoeff(len_segs(0, 0) +1, len_segs(1, 0) + Slacks_vec + 1, Slacks_vec, "index");
    }
    Int_Coeff[nSeg[0]-1] = C_InterpCoeff(len_segs(0, nSeg[0]-1) +1, len_segs(1, nSeg[0]-1) + Slacks_vec + 1, Slacks_vec, "coeff");
    Int_Index[nSeg[0]-1] = C_InterpCoeff(len_segs(0, nSeg[0]-1) +1, len_segs(1, nSeg[0]-1) + Slacks_vec + 1, Slacks_vec, "index");
  } else {
    // seg length > 1
    for(int i=0; i < nSeg[0]; ++i){
      Int_Coeff[i] = C_InterpCoeff(len_segs(0, i) +1, len_segs(1, i) + Slacks_vec + 1, Slacks_vec, "coeff");
      Int_Index[i] = C_InterpCoeff(len_segs(0, i) +1, len_segs(1, i) + Slacks_vec + 1, Slacks_vec, "index");
    }
  }
  
  //   #### Dynamic Programming Section ####
  
  NumericVector table_index(Bounds.cols()+1);
  table_index[0] = 0;
  //  table_index = cumsum(c(0, Bounds[2, ] - Bounds[1, ] + 1) )
  for(int i=1; i < table_index.size(); ++i){
    table_index[i] = Bounds(1, i-1) - Bounds(0, i-1) + 1;
    table_index[i] = table_index[i] + table_index[i - 1];
  }
  NumericMatrix Table(3, table_index[nSeg[0] + 1]);
  for(int i=1; i < (Table.cols()); ++i){
    Table(1,i) = R_NegInf; // from https://teuder.github.io/rcpp4everyone_en/240_na_nan_inf.html
  }
  
  // below replaces below R; 
  //    seems to work but I do not completely understand what it is doing 
  //    or if original works how it was meant to...
  // for (i in seq(1, nSeg + 1)) {
  //  v <- seq(Bounds[1, i], Bounds[2, i])
  //  Table[1, seq(table_index[i] +1, table_index[i +1])  ] <-v
  //}
  for (int i=0; i < table_index.length()-1; ++i) {
      for(int j = 0; j < (Bounds(1, i) - Bounds(0, i))+1; ++j){
        Table(0, table_index[i]+j) = Bounds(0, i)+j;
      }
  }
  
  // Forward phase
  // ignoring R* time1 <- Sys.time()
  
  //next bit is the big bit....
  
  //output 
  //passing to list so I do not have keep resetting output declaration in top line...
  //be a list at end anyway....
  
  // test start
  // List test(nSeg[0]);
  // int tmp =0;
  
  for (int i = 0; i<nSeg[0]; ++i){
    IntegerVector a = Slacks_vec + len_segs(1, i);
    int b = table_index[i] + 1 - Bounds(0, i);
    int c = len_segs(0, i) + 1;
    int counting = 1;
    int node_z = table_index(i + 2);
    int node_a = table_index(i + 1) + 1;
    NumericMatrix bound_k_table(2, node_z - node_a + 1); // these are not NAs...
    // this should be simpler ???
    // R * int_index_seg <- t(Int_Index[[i]]) - (len_segs[2, i] +1)
    NumericMatrix int_index_tmp = Int_Index[i];
    NumericMatrix int_index_seg(int_index_tmp.cols(), int_index_tmp.rows());
    for(int ii =0; ii < int_index_tmp.cols(); ++ii){
      for(int jj =0; jj < int_index_tmp.rows(); ++jj){
        int_index_seg(ii,jj) = int_index_tmp(jj,ii);
      }
    }
    int_index_seg = int_index_seg - (len_segs(1,i) + 1);
    // this above... ???
    // R * int_coeff_seg <- t(Int_Coeff[[i]])
    NumericMatrix int_coeff_tmp = Int_Coeff[i];
    NumericMatrix int_coeff_seg(int_coeff_tmp.cols(), int_coeff_tmp.rows());
    for(int ii =0; ii < int_coeff_tmp.cols(); ++ii){
      for(int jj =0; jj < int_coeff_tmp.rows(); ++jj){
        int_coeff_seg(ii,jj) = int_coeff_tmp(jj,ii);
      }
    }
    NumericVector TSeg(bT[i + 1]-bT[i]+1); 
    for(int ii =0; ii < TSeg.length(); ++ii){
      TSeg[ii] = Ta[bT[i]-1+ii];
    }
    NumericVector TSeg_centered = TSeg - (sum(TSeg)/TSeg.length());
    // from https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/norm
    // R* Norm_TSeg_cen  <- norm(TSeg_centered, type = '2')
    NumericVector Norm_TSeg_cen(1); 
    Norm_TSeg_cen = sqrt(sum(TSeg_centered* TSeg_centered));
    for (int j =node_a; j < node_z+1; ++j)  {
      IntegerVector prec_nodes = Table(0, j-1) - a; 
      LogicalVector allowed_arcs = ( prec_nodes >= Bounds(0, i) ) & ( prec_nodes <= Bounds(1, i) );
      IntegerVector nodes_tablePointer = prec_nodes[allowed_arcs];
      nodes_tablePointer = nodes_tablePointer + b; // MIGHT BE 1 OUT because this R-to-C...
      int n_aa = sum(allowed_arcs);
      
      if(n_aa  != 0 ){
        // got to line 206...
        // Index_Node <- Table[1, j] + int_index_seg[, allowed_arcs]
        NumericMatrix Index_Node(int_index_seg.rows(), n_aa);
        
        int kk = 0;
        for(int ii=0; ii<allowed_arcs.length(); ++ii){
          if(allowed_arcs[ii]){
            for(int jj=0; jj < int_index_seg.rows(); ++jj){
              Index_Node(jj, kk) = int_index_seg(jj, ii);
            }
            kk = kk + 1;
          }
        }
        Index_Node = Index_Node + Table(0, j-1);
        
        // coeff_b <- sapply(int_coeff_seg[, allowed_arcs], function(x) x)
        NumericMatrix Index_tmp(int_index_seg.rows(), n_aa);
        kk = 0;
        for(int ii=0; ii<allowed_arcs.length(); ++ii){
          if(allowed_arcs[ii]){
            for(int jj=0; jj < int_coeff_seg.rows(); ++jj){
              Index_tmp(jj, kk) = int_coeff_seg(jj, ii);
            }
            kk = kk + 1;
          }
        }
        NumericVector coeff_b(Index_tmp.begin(), Index_tmp.end());
        
        // Xi_seg <- X[Index_Node]
        NumericMatrix Xi_seg_sp(Index_Node.rows(), Index_Node.cols()); 
        for(int ii=0; ii<Index_Node.ncol(); ++ii){
          for(int jj=0; jj<Index_Node.nrow(); ++jj){
            Xi_seg_sp(jj, ii) = X(0, Index_Node(jj, ii)-1); // NB this only does first row of matrix
          }
        }
        // Xi_diff <- Xdiff[Index_Node]
        NumericMatrix Xi_diff(Index_Node.rows(), Index_Node.cols()); 
        for(int ii=0; ii<Index_Node.ncol(); ++ii){
          for(int jj=0; jj<Index_Node.nrow(); ++jj){
            Xi_diff(jj, ii) = Xdiff(0, Index_Node(jj, ii)-1);  // NB this only does first row of matrix
          }
        }
        
        // toreshape <- matrix(c(coeff_b, Xi_diff), nrow = 2, byrow = T)
        // toreshape <- apply(toreshape, 2, prod)  + Xi_seg
        NumericVector reshape1(Xi_diff.begin(), Xi_diff.end());
        for(int ii=0; ii<coeff_b.length(); ++ii){
            reshape1[ii] = coeff_b[ii] * reshape1[ii];
        }
        reshape1 = reshape1 + Xi_seg_sp;
        // Xi_seg <- matrix(toreshape, nrow = c, ncol = n_aa*dimX[1] )
        NumericMatrix  Xi_seg(c, n_aa * dimX[0], reshape1.begin() );
        // Xi_Seg_mean <- colSums(Xi_seg)/nrow(Xi_seg)
        // Norm_Xi_Seg_cen <- sqrt(colSums(Xi_seg^2) - nrow(Xi_seg) * Xi_Seg_mean^2)
        NumericVector Xi_Seg_mean(Xi_seg.cols());
        for(int ii=0; ii<Xi_seg.cols(); ++ii){
          Xi_Seg_mean[ii] = mean(Xi_seg.column(ii));
        }
        NumericVector Norm_Xi_Seg_cen(Xi_seg.cols());
        for(int ii=0; ii<Norm_Xi_Seg_cen.length(); ++ii){
          Norm_Xi_Seg_cen[ii] = sum(Xi_seg.column(ii) * Xi_seg.column(ii) );
        }
        Norm_Xi_Seg_cen = sqrt( Norm_Xi_Seg_cen - Xi_seg.rows() * (Xi_Seg_mean * Xi_Seg_mean) );  
        
        //  CCs_Node <- ( as.numeric(TSeg_centered) %*% Xi_seg ) /
        //  (Norm_TSeg_cen %*% Norm_Xi_Seg_cen)
        // NumericMatrix CCs_Node = as<NumericMatrix>(as<NumericVector>(TSeg_centered)) * as<NumericMatrix>(Xi_seg); //  / (Norm_TSeg_cen * Norm_Xi_Seg_cen) ;  
        NumericMatrix CCs_Node(1, Xi_seg.cols());
        for(int ii=0; ii < Xi_seg.cols(); ++ii){
          CCs_Node(0,ii) = sum( as<NumericVector>(TSeg_centered) * Xi_seg.column(ii) );
        }
        // (Norm_TSeg_cen * Norm_Xi_Seg_cen)
        // CCs_Node <- ifelse(is.finite(CCs_Node), CCs_Node, 0)
        NumericVector CCs_Node_b = Norm_TSeg_cen[0] * Norm_Xi_Seg_cen;
        for(int ii=0; ii < CCs_Node.ncol(); ++ii){
          if(std::isfinite((CCs_Node(0,ii) / CCs_Node_b[ii]))){
            CCs_Node(0,ii) = (CCs_Node(0,ii) / CCs_Node_b[ii]);
          } else {
            CCs_Node(0,ii) = 0;
          }
        }
        //if(Options[[2]]){
        //  Cost_Fun = matrix(Table[2, nodes_tablePointer], nrow = n_aa ) + CCs_Node
        //} else {
        //  Cost_Fun <- matrix(Table[2, nodes_tablePointer], nrow = n_aa) + CCs_Node^Options[[2]]
        //}
        CCs_Node.attr("dim") = Dimension(n_aa, dimX[0]);
        NumericMatrix Cost_Fun(nodes_tablePointer.length(), 1);
        for(int ii=0; ii<nodes_tablePointer.length(); ++ii){
          if(Options[1]==1){
            Cost_Fun(ii,0) = Table(1, nodes_tablePointer[ii]-1) + CCs_Node(0,ii);
          } else {
            Cost_Fun(ii,0) = Table(1, nodes_tablePointer[ii]-1) + pow(CCs_Node(0,ii),Options[1]);
          }
        }
        
        //ind <- max(Cost_Fun)
        //pos <- match(ind, Cost_Fun)
        double ind = max(Cost_Fun);
        double pos = which_max(Cost_Fun);
        bound_k_table(0, counting-1) = ind;
        bound_k_table(1, counting-1) = nodes_tablePointer[pos];  // remember 0 not NAs...
        
        counting = counting + 1;

        // test(i) = bound_k_table;
        // test end
         
      }
      // the counter for testing within for loops  
      //tmp=tmp+1;
    } 
    //Table[2:3, seq(node_a, node_z)] <- bound_k_table
    for(int ii=0; ii< bound_k_table.cols(); ++ii){
      Table(1, node_a+ii-1) = bound_k_table(0, ii);
      Table(2, node_a+ii-1) = bound_k_table(1, ii);
    }
    
  }
  
  // return test;
  
  // time2 <- Sys.time()
  // elapsed <- time2 - time1
  // ignoring time keeping 
  
  //for (i in seq(1, dimX[1])) {
  //  Pointer <- ncol(Table)
  //  Warping[i, nSeg + 1] <- dimX[2]
  //  for (j in seq(nSeg, 1, -1)) {
  //    Pointer <- Table[3, Pointer]
  //    Warping[i, j] <- Table[1, Pointer]
  //  }
  //}
  
  for(int i=0; i < dimX[0]; ++i){
    int Pointer = Table.cols() -1;
    Warping(i, nSeg[0]) = dimX[1]; 
    for(int j= nSeg[0]-1; j>-1; --j){
      Pointer = Table(2, Pointer)-1; 
      Warping(i,j) = Table(0, Pointer);
    }
  }
  
  //Xwarped <- NULL
  //for (i in seq(1, nSeg)) {
  //  indT <- seq(bT[i], bT[i + 1])
  //  lenT <- bT[i + 1] - bT[i]
  //  for (j in seq(1, dimX[1])) {
  //    indX <- seq(Warping[j, i], Warping[j, i + 1])
  //    lenX <- Warping[j, i + 1] - Warping[j, i]
  //    Xwarped[indT] <- approx(x = indX - Warping[j, i] + 1,
  //                            y = X[indX],
  //                                 xout = seq(0, lenT)/lenT * lenX + 1 )$y
  //  }
  //}
  
  Function R_approx( "approx" );
  // test(0) =  Xwarped; // R_approx(Warping, bT, seq(1, X.length()));
  
  //testing
  //List out(5);
  
  for(int i=0; i<bT.length()-1; ++i){
    double lenT = bT[i+1] - bT[i];
    NumericVector indT(lenT+1); 
    for(int ii=0; ii<lenT+1; ++ii){
      indT[ii] = bT[i] + ii;
    }
    for(int j=0; j<dimX[0]; ++j){
      double lenX = Warping(j, i + 1) - Warping(j, i);
      NumericVector indX(lenX+1); 
      for(int ii=0; ii<lenX+1; ++ii){
        indX[ii] = Warping(j, i) + ii;
      }
      NumericVector y_temp(indX.length());
      for(int ii=0; ii<y_temp.length(); ++ii){
        y_temp[ii] = X(j, indX[ii]-1);
      }
      NumericVector new_temp(lenT+1);
      for(int ii=0; ii<new_temp.length(); ++ii){
        new_temp[ii] = ii; 
      }
      new_temp = new_temp / (lenT) * (lenX) + 1;
      //out(i) = new_temp;
      List app_out = R_approx(indX - Warping(j, i) + 1, y_temp, new_temp );
      NumericVector out_temp = app_out["y"];
      for(int ii=0; ii<indT.length(); ++ii){
        Xwarped(j, indT(ii)-1) = out_temp[ii]; 
      }
    }
  }

  //return(list(bT = bT, Warping = Warping, nSeg=nSeg,
  //            Xwarped = Xwarped))
  List out = List::create(Named("bT", bT),
                          Named("Warping", Warping),
                          Named("nSeg", nSeg),
                          Named("Xwarped"), Xwarped);
  return out;

}

// C version of cow... in R  
// could improve/tidy this a lot but it works!!!!
// hurts my head thinking about this...

// this this should generate following
// /*** R
// > cow(1:10, 1:10, Seg=3, Slack=1, Options = c(2,1,0,1))
// > C_cow(1:10, matrix(1:10), Seg=3, Slack=1, Options = c(2,1,0,1)) 
// [1] "Segments: 3 Points: 3"
// $bT
// [1]  1  3  5  7 10
// 
// $Warping
// [,1] [,2] [,3] [,4] [,5]
// [1,]    1    4    7    8   10
// 
// $nSeg
// [1] 4
// 
// $Xwarped
// [1]  1.000000  2.500000  4.000000  5.500000  7.000000  7.500000  8.000000  8.666667  9.333333 10.000000
// */

// test 
// Rcpp::compileAttributes(); devtools::load_all(); C_cow(1:10, matrix(1:10, nrow=1), Seg=3, Slack=1, Options = c(1,1,0,1)); cow(1:10, 1:10, Seg=3, Slack=1, Options = c(1,1,0,1))
