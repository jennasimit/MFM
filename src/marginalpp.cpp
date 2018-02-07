#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
int stroverlap(const IntegerVector& x,
	       const IntegerVector& y) {
  int nx=x.size();
  int ny=y.size();
  if(nx==0)
    return 0;
  if(ny==0)
    return 0;
  for(int ix=0; ix<nx; ix++)
    for(int iy=0; iy<ny; iy++)
      if(x[ix]==y[iy])
	return 1;
  return 0;
}

// returns x/sum(x)
// Rccp::export
NumericVector sum1(const NumericVector& x) {
  const int n=x.size();
  double sum=0.0;
  NumericVector y(n);
  for(int i=0; i<n; i++)
    sum+=x(i);
  for(int i=0; i<n; i++)
    y(i)=x(i)/sum;
  return(y);
}


// [[Rcpp::export]]
List calcQpair(const List S1, // model matrix - columns are models, rows are SNPs
	       const List S2, // model matrix - columns are models, rows are SNPs
	       const NumericVector& pp1,
	       const NumericVector& pp2) { // pp for each model for disease 2
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  // const int nsnp = M.nrow();
  NumericVector Q1(nmod1);
  NumericVector Q2(nmod2);
  for(int i=0; i<nmod1; i++) {
    for(int j=0; j<nmod2; j++) {
      int idx=stroverlap( S1[i], S2[j]);
      if(idx > 0) {
	Q1(i) = Q1(i) + pp2(j);
	Q2(j) = Q2(j) + pp1(i);
      }
    }
  }

  // NB make each Q sum to 1
  List Q = List::create(sum1(Q1), sum1(Q2));
  return(Q);
}

NumericMatrix bind2(const NumericVector A,
		    const NumericVector B) {
  return(Rcpp::cbind(A,B));
}
  
NumericMatrix bind3(const NumericVector A,
		    const NumericVector B,
		    const NumericVector C) {
  return(Rcpp::cbind(A,B,C));
}
  
NumericMatrix bind4(const NumericVector A,
		    const NumericVector B,
		    const NumericVector C,
		    const NumericVector D) {
  return(Rcpp::cbind(A,B,C,D));
}
  
// [[Rcpp::export]]
List newcalcQ2(const List S1,
	    const List S2,
	    const NumericVector& pp1,
	    const NumericVector& pp2) { // pp for each model for disease 2
  List Q12 = calcQpair(S1,S2,pp1,pp2);
  Q12.names() = CharacterVector::create("1", "2");
  return(wrap(Q12));
}

// [[Rcpp::export]]
List newcalcQ3(const List S1,
	    const List S2,
	    const List S3, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3) { // pp for each model for disease 2
  List Q12 = calcQpair(S1,S2,pp1,pp2);
  List Q23 = calcQpair(S2,S3,pp2,pp3);
  List Q13 = calcQpair(S1,S3,pp1,pp3);
    
  // const int nmod1 = pp1.size();
  // const int nmod2 = pp2.size();
  // const int nmod3 = pp3.size();
  NumericMatrix Q1 = bind2(Q12(0),Q13(0)); // Q for trait 1 | 2
  NumericMatrix Q2 = bind2(Q12(1),Q23(0));
  NumericMatrix Q3 = bind2(Q13(1),Q23(1));

  List Q = List::create(_["1"] = Q1,
			_["2"] = Q2,
			_["3"] = Q3);
  return(wrap(Q));
}

// [[Rcpp::export]]
List newcalcQ4(const List S1,
	    const List S2,
	    const List S3, // model matrix - columns are models, rows are SNPs
	    const List S4, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3,
	    const NumericVector& pp4) { // pp for each model for disease 2
  List Q12 = calcQpair(S1,S2,pp1,pp2);
  List Q13 = calcQpair(S1,S3,pp1,pp3);
  List Q14 = calcQpair(S1,S4,pp1,pp4);
  List Q23 = calcQpair(S2,S3,pp2,pp3);
  List Q24 = calcQpair(S2,S4,pp2,pp4);
  List Q34 = calcQpair(S3,S4,pp3,pp4);
    
  // const int nmod1 = pp1.size();
  // const int nmod2 = pp2.size();
  // const int nmod3 = pp3.size();
  // const int nmod4 = pp3.size();
  NumericMatrix Q1 = bind3(Q12[0],Q13[0],Q14[0]); // Q for trait 1 | 2
  NumericMatrix Q2 = bind3(Q12(1),Q23(0),Q24(0));
  NumericMatrix Q3 = bind3(Q13(1),Q23(1),Q34(0));
  NumericMatrix Q4 = bind3(Q14[1],Q24[1],Q34[1]);

  List Q = List::create(_["1"] = Q1,
			_["2"] = Q2,
			_["3"] = Q3,
			_["4"] = Q4);

  return(wrap(Q));
}


