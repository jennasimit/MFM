#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector stroverlap(const IntegerVector& x,
	       const IntegerVector& y) {
  int nx=x.size();
  int ny=y.size();
  IntegerVector ret = IntegerVector::create(0,nx,ny);
  if(nx==0)
    return ret;
  if(ny==0)
    return ret;
  for(int ix=0; ix<nx; ix++)
    for(int iy=0; iy<ny; iy++)
      if(x[ix]==y[iy]) {
	ret[0]=1;
	return(ret);
      }
  return ret;
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
	       const NumericVector& pp2, // pp for each model for disease 2
	       const NumericMatrix& tau,
	       const double kappa) { 
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  // const int ntau = tau.cols();
  // const int nsnp = M.nrow();
  NumericVector Q1(nmod1);
  NumericVector Q2(nmod2);
  double k=0.0;
  for(int i=0; i<nmod1; i++) {
    for(int j=0; j<nmod2; j++) {
      if(i==0 | j==0) { // at least one is null model -> tau=1 & int is empty
	k=1.0;
      } else {
	IntegerVector idx=stroverlap( S1[i], S2[j]);
	int n1=idx(1);
	int n2=idx(2);
	if(idx(0) == 0) {
	  k = tau(n1,n2);
	} else {
	  k = kappa * tau(n1,n2);
	}
      }
      Q1(i) = Q1(i) + pp2(j) * k;
      Q2(j) = Q2(j) + pp1(i) * k;
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
NumericMatrix bind5(const NumericVector A,
		    const NumericVector B,
		    const NumericVector C,
		    const NumericVector D,
		    const NumericVector E) {
  return(Rcpp::cbind(A,B,C,D,E));
}
   
// [[Rcpp::export]]
List newcalcQ2(const List S1,
	    const List S2,
	    const NumericVector& pp1,
	       const NumericVector& pp2,
	       const NumericMatrix& tau,
	       const double kappa) { // pp for each model for disease 2
  List Q12 = calcQpair(S1,S2,pp1,pp2,tau,kappa);
  Q12.names() = CharacterVector::create("1", "2");
  return(wrap(Q12));
}


// [[Rcpp::export]]
List newcalcQ3(const List S1,
	    const List S2,
	    const List S3, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3,
	       const NumericMatrix& tau,
	       const double kappa) { // pp for each model for disease 2
  List Q12 = calcQpair(S1,S2,pp1,pp2,tau,kappa);
  List Q23 = calcQpair(S2,S3,pp2,pp3,tau,kappa);
  List Q13 = calcQpair(S1,S3,pp1,pp3,tau,kappa);
    
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
	    const NumericVector& pp4,
	       const NumericMatrix& tau,
	       const double kappa) { // pp for each model for disease 2
  List Q12 = calcQpair(S1,S2,pp1,pp2,tau,kappa);
  List Q13 = calcQpair(S1,S3,pp1,pp3,tau,kappa);
  List Q14 = calcQpair(S1,S4,pp1,pp4,tau,kappa);
  List Q23 = calcQpair(S2,S3,pp2,pp3,tau,kappa);
  List Q24 = calcQpair(S2,S4,pp2,pp4,tau,kappa);
  List Q34 = calcQpair(S3,S4,pp3,pp4,tau,kappa);
    
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


// [[Rcpp::export]]
List newcalcQ5(const List S1,
	    const List S2,
	    const List S3, // model matrix - columns are models, rows are SNPs
	    const List S4, // model matrix - columns are models, rows are SNPs
	    const List S5, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3,
	    const NumericVector& pp4,
	    const NumericVector& pp5,
	       const NumericMatrix& tau,
	       const double kappa) { // pp for each model for disease 2
  List Q12 = calcQpair(S1,S2,pp1,pp2,tau,kappa);
  List Q13 = calcQpair(S1,S3,pp1,pp3,tau,kappa);
  List Q14 = calcQpair(S1,S4,pp1,pp4,tau,kappa);
  List Q15 = calcQpair(S1,S5,pp1,pp5,tau,kappa);
  List Q23 = calcQpair(S2,S3,pp2,pp3,tau,kappa);
  List Q24 = calcQpair(S2,S4,pp2,pp4,tau,kappa);
  List Q25 = calcQpair(S2,S5,pp2,pp5,tau,kappa);
  List Q34 = calcQpair(S3,S4,pp3,pp4,tau,kappa);
  List Q35 = calcQpair(S3,S5,pp3,pp5,tau,kappa);
  List Q45 = calcQpair(S4,S5,pp4,pp5,tau,kappa);
    
  // const int nmod1 = pp1.size();
  // const int nmod2 = pp2.size();
  // const int nmod3 = pp3.size();
  // const int nmod4 = pp3.size();
  NumericMatrix Q1 = bind4(Q12[0],Q13[0],Q14[0],Q15[0]); // Q for trait 1 | 2
  NumericMatrix Q2 = bind4(Q12(1),Q23(0),Q24(0),Q25[0]);
  NumericMatrix Q3 = bind4(Q13(1),Q23(1),Q34(0),Q35[0]);
  NumericMatrix Q4 = bind4(Q14[1],Q24[1],Q34[1],Q45[0]);
  NumericMatrix Q5 = bind4(Q15[1],Q25[1],Q35[1],Q45[1]);

  List Q = List::create(_["1"] = Q1,
			_["2"] = Q2,
			_["3"] = Q3,
			_["4"] = Q4,
			_["5"] = Q5);

  return(wrap(Q));
}




// [[Rcpp::export]]
List newcalcQ6(const List S1,
	    const List S2,
	    const List S3, // model matrix - columns are models, rows are SNPs
	    const List S4, // model matrix - columns are models, rows are SNPs
	    const List S5, // model matrix - columns are models, rows are SNPs
	    const List S6, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3,
	    const NumericVector& pp4,
	    const NumericVector& pp5,
	    const NumericVector& pp6,
	       const NumericMatrix& tau,
	       const double kappa) { // pp for each model for disease 2
  List Q12 = calcQpair(S1,S2,pp1,pp2,tau,kappa);
  List Q13 = calcQpair(S1,S3,pp1,pp3,tau,kappa);
  List Q14 = calcQpair(S1,S4,pp1,pp4,tau,kappa);
  List Q15 = calcQpair(S1,S5,pp1,pp5,tau,kappa);
  List Q16 = calcQpair(S1,S6,pp1,pp6,tau,kappa);
  List Q23 = calcQpair(S2,S3,pp2,pp3,tau,kappa);
  List Q24 = calcQpair(S2,S4,pp2,pp4,tau,kappa);
  List Q25 = calcQpair(S2,S5,pp2,pp5,tau,kappa);
  List Q26 = calcQpair(S2,S6,pp2,pp6,tau,kappa);
  List Q34 = calcQpair(S3,S4,pp3,pp4,tau,kappa);
  List Q35 = calcQpair(S3,S5,pp3,pp5,tau,kappa);
  List Q36 = calcQpair(S3,S6,pp3,pp6,tau,kappa);
  List Q45 = calcQpair(S4,S5,pp4,pp5,tau,kappa);
  List Q46 = calcQpair(S4,S6,pp4,pp6,tau,kappa);
  List Q56 = calcQpair(S5,S6,pp5,pp6,tau,kappa);
    
  // const int nmod1 = pp1.size();
  // const int nmod2 = pp2.size();
  // const int nmod3 = pp3.size();
  // const int nmod4 = pp3.size();
  NumericMatrix Q1 = bind5(Q12[0],Q13[0],Q14[0],Q15[0],Q16[0]); // Q for trait 1 | 2
  NumericMatrix Q2 = bind5(Q12(1),Q23(0),Q24(0),Q25[0],Q26[0]);
  NumericMatrix Q3 = bind5(Q13(1),Q23(1),Q34(0),Q35[0],Q36[0]);
  NumericMatrix Q4 = bind5(Q14[1],Q24[1],Q34[1],Q45[0],Q46[0]);
  NumericMatrix Q5 = bind5(Q15[1],Q25[1],Q35[1],Q45[1],Q56[0]);
  NumericMatrix Q6 = bind5(Q16[1],Q26[1],Q36[1],Q46[1],Q56[1]);

  List Q = List::create(_["1"] = Q1,
			_["2"] = Q2,
			_["3"] = Q3,
			_["4"] = Q4,
			_["5"] = Q5,
			_["6"] = Q6);

  return(wrap(Q));
}

