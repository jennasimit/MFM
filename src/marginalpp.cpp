#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
int modoverlap(const arma::ivec& x,
	       const arma::ivec& y) {
  int idx = sum(x % y);
  if(idx > 1)
    idx = 1;
  return idx;
}

// [[Rcpp::export]]
List calcQ2(const arma::imat& M1, // model matrix - columns are models, rows are SNPs
	    const arma::imat& M2, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2) { // pp for each model for disease 2
  Rcpp::Rcout << "size of M: " << size(M1) << std::endl;
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  // const int nsnp = M.nrow();
  NumericVector Q1(nmod1);
  NumericVector Q2(nmod2);
  for(int i=0; i<nmod1; i++) {
    for(int j=0; j<nmod2; j++) {
      int idx=modoverlap( M1.col(i), M2.col(j));
      if(idx > 0) {
	Q1(i) = Q1(i) + pp2(j);
	Q2(j) = Q2(j) + pp1(i);
      }
    }
  }
  List Q(2);
  Q[0] = Q1;
  Q[1] = Q2;
  return(wrap(Q));
}


// [[Rcpp::export]]
List calcQ3(const arma::imat& M1,
	    const arma::imat& M2,
	    const arma::imat& M3, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3) { // pp for each model for disease 2
  Rcpp::Rcout << "size of M: " << size(M1) << std::endl;
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  const int nmod3 = pp3.size();
  NumericMatrix Q1(nmod1,3);
  NumericMatrix Q2(nmod2,3);
  NumericMatrix Q3(nmod3,3);
  for(int i=0; i<nmod1; i++) {
    for(int j=0; j<nmod2; j++) {
      int idx1=sum( M1.col(i) % M2.col(j));
      for(int k=0; k<nmod3; k++) {
	Rcpp::Rcout << i << ' ' << j << ' ' << k << std::endl;
	int idx=idx1 +
	  modoverlap(M1.col(i), M3.col(k)) +
	  modoverlap(M2.col(j), M3.col(k));
	if(idx > 0) {
	  Q1(i,idx) = Q1(i,idx) + pp2(j) * pp3(k);
	  Q2(j,idx) = Q2(j,idx) + pp1(i) * pp3(k);
	  Q3(j,idx) = Q3(j,idx) + pp1(i) * pp2(k);
	}
      }
    }
  }
  List Q(3);
  Q[0] = Q1;
  Q[1] = Q2;
  Q[2] = Q3;
  return(wrap(Q));
}

  
// List calcQ(List M, // list of model matrices, each will be a const arma::imat& - columns are models, rows are SNPs
// 		    List pp) { // List of pp for each model for each disease
//   const int nM = M.size();
//   const int npp = pp.size();
//   Rcpp::Rcout << "number of things: " << nM << " or " << npp << std::endl;
//   List Q(nM);
//   // // const int nsnp = M.nrow();
//   // NumericVector S(2);
//   // NumericVector Q(nmod);
//   for(int im = 0; im<nM; im++) {
//     Q[im] = arma::zeros(M[im].size()); 
//   // //   S(i) = 0.0;
  
//   // for(int i =0; i<nmod; i++) {
//   //   for(int j=0; j<nmod; j++) {
//   //     int idx=sum( M.col(i) % M.col(j));
//   //     if(idx > 1)
//   // 	Q(i) = Q(j) + pp2(j);
//   //   }
//   // }
//   }
//   return(Q);
// }


