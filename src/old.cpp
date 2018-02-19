#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/*
// [[Rcpp::export]]
int modoverlap(const arma::ivec& x,
	       const arma::ivec& y) {
  int idx = sum(x % y);
  if(idx > 1)
    idx = 1;
  return idx;
}


// slower
// [[Rcpp::export]]
int strint(const IntegerVector& x,
	   const IntegerVector& y) {
  IntegerVector insect = intersect(x,y);
  if(insect.size()>0)
    return 1;
  return 0;
}

// [[Rcpp::export]]
List calcQ2(const List S1, // model matrix - columns are models, rows are SNPs
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
  List Q = List::create(_["1"] = Q1, _["2"] = Q2);
  return(wrap(Q));
}



// [[Rcpp::export]]
List calcQ2_models(const arma::imat& M1, // model matrix - columns are models, rows are SNPs
	    const arma::imat& M2, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2) { // pp for each model for disease 2
  // Rcpp::Rcout << "size of M: " << size(M1) << std::endl;
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
List calcQ3(const List S1,
	    const List S2,
	    const List S3, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3) { // pp for each model for disease 2
  // Rcpp::Rcout << "size of M: " << size(M1) << std::endl;
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  const int nmod3 = pp3.size();
  NumericMatrix Q1(nmod1,3);
  NumericMatrix Q2(nmod2,3);
  NumericMatrix Q3(nmod3,3);
  for(int i1=0; i1<nmod1; i1++) {
     for(int i2=0; i2<nmod2; i2++) {
       int idx1=stroverlap(S1[i1], S2[i2]) - 1;
       double pp12 = pp1(i1) * pp2(i2);
       for(int i3=0; i3<nmod3; i3++) {
   	int idx=idx1 +
   	  stroverlap(S1[i1], S3[i3]) +
   	  stroverlap(S2[i2], S3[i3]);
	if(idx > 0) {
	  // idx = idx - 1;
	  Q1(i1,idx) = Q1(i1,idx) + pp2(i2) * pp3(i3);
  	  Q2(i2,idx) = Q2(i2,idx) + pp1(i1) * pp3(i3);
  	  Q3(i3,idx) = Q3(i3,idx) + pp12;
  	}
       }
     }
  }
  List Q;
  Q["1"] = Q1;
  Q["2"] = Q2;
  Q["3"] = Q3;
  return(Q);
}



// [[Rcpp::export]]
List intcalcQ3(const List S1,
	    const List S2,
	    const List S3, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3) { // pp for each model for disease 2
  List Q12 = calcQpair(S1,S2,pp1,pp2);
  List Q23 = calcQpair(S2,S3,pp2,pp3);
  List Q13 = calcQpair(S1,S3,pp1,pp3);
    
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  const int nmod3 = pp3.size();
  // const int nsnp = M.nrow();
  NumericVector Q1 = Q12[0];
  NumericVector Q2 = Q12[1];
  NumericVector Q3 = Q13[1];

  NumericVector tmp = Q13[0];
  for(int i=0; i<nmod1; i++)
    Q1(i)+= tmp(i);
  tmp = Q23[0];
  for(int i=0; i<nmod2; i++)
    Q2(i)+= tmp(i);
  tmp = Q23[1];
  for(int i=0; i<nmod3; i++)
    Q3(i)+= tmp(i);

  List Q(3);
  Q["1"] = Q1;
  Q["2"] = Q2;
  Q["3"] = Q3;
  return(wrap(Q));
}



// [[Rcpp::export]]
List intcalcQ4(const List S1,
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
    
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  const int nmod3 = pp3.size();
  const int nmod4 = pp4.size();
  // const int nsnp = M.nrow();
  NumericVector Q1 = Q12[0];
  NumericVector Q2 = Q12[1];
  NumericVector Q3 = Q13[1];
  NumericVector Q4 = Q14[1];

  NumericVector a,b;
  a = Q13[0];
  b = Q14[0];
  for(int i=0; i<nmod1; i++)
    Q1(i)+= a(i) + b(i);
  a = Q23[0];
  b = Q24[0];
  for(int i=0; i<nmod2; i++)
    Q2(i)+= a(i) + b(i);
  a = Q23[1];
  b = Q34[0];
  for(int i=0; i<nmod3; i++)
    Q3(i)+= a(i) + b(i);
  a = Q24[1];
  b = Q34[1];
  for(int i=0; i<nmod4; i++)
    Q4(i)+= a(i) + b(i);

  List Q(4);
  Q["1"] = Q1;
  Q["2"] = Q2;
  Q["3"] = Q3;
  Q["4"] = Q4;
  return(wrap(Q));
}



// [[Rcpp::export]]
List calcQ3log(const List S1,
	    const List S2,
	    const List S3, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3) { // pp for each model for disease 2
  // Rcpp::Rcout << "size of M: " << size(M1) << std::endl;
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  const int nmod3 = pp3.size();
  NumericMatrix Q1(nmod1,3);
  NumericMatrix Q2(nmod2,3);
  NumericMatrix Q3(nmod3,3);
  for(int i1=0; i1<nmod1; i1++) {
     for(int i2=0; i2<nmod2; i2++) {
       int idx1=stroverlap(S1[i1], S2[i2]) - 1;
       double ep12 = exp(pp1(i1) + pp2(i2));
       for(int i3=0; i3<nmod3; i3++) {
   	int idx=idx1 +
   	  stroverlap(S1[i1], S3[i3]) +
   	  stroverlap(S2[i2], S3[i3]);
	if(idx > 0) {
	  // idx = idx - 1;
	  Q1(i1,idx) = Q1(i1,idx) + exp(pp2(i2) + pp3(i3));
  	  Q2(i2,idx) = Q2(i2,idx) + exp(pp1(i1) + pp3(i3));
  	  Q3(i3,idx) = Q3(i3,idx) + ep12;
  	}
       }
     }
  }
  List Q;
  Q["1"] = Q1;
  Q["2"] = Q2;
  Q["3"] = Q3;
  return(Q);
}



// [[Rcpp::export]]
List calcQ4(const List S1,
	    const List S2,
	    const List S3, // model matrix - columns are models, rows are SNPs
	    const List S4, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3,
	    const NumericVector& pp4) { // pp for each model for disease 2
  // Rcpp::Rcout << "size of M: " << size(M1) << std::endl;
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  const int nmod3 = pp3.size();
  const int nmod4 = pp4.size();
  NumericMatrix Q1(nmod1,6);
  NumericMatrix Q2(nmod2,6);
  NumericMatrix Q3(nmod3,6);
  NumericMatrix Q4(nmod4,6);
  for(int i1=0; i1<nmod1; i1++) {
     for(int i2=0; i2<nmod2; i2++) {
       int idx1=stroverlap(S1[i1], S2[i2]) - 1;
       double pp12 = pp1(i1) * pp2(i2);
       for(int i3=0; i3<nmod3; i3++) {
	 int idx2=idx1 +
	   stroverlap(S1[i1], S3[i3]) +
	   stroverlap(S2[i2], S3[i3]);
	 double pp23 = pp2(i2) * pp3(i3);
	 double pp13 = pp1(i1) * pp3(i3);
	 double pp123 = pp1(i1) * pp2(i2) * pp3(i3);
	 for(int i4=0; i4<nmod4; i4++) {
	   int idx=idx2 +
	     stroverlap(S1[i1], S4[i4]) +
	     stroverlap(S2[i2], S4[i4]) +
	     stroverlap(S3[i3], S4[i4]);
	   if(idx > 0) {
	     // idx = idx - 1;
	     Q1(i1,idx) = Q1(i1,idx) + pp23 * pp4(i4);
	     Q2(i2,idx) = Q2(i2,idx) + pp13 * pp4(i4);
	     Q3(i3,idx) = Q3(i3,idx) + pp12 * pp4(i4);
	     Q4(i4,idx) = Q4(i4,idx) + pp123;
	   }
	 }
       }
     }
  }
  List Q;
  Q["1"] = Q1;
  Q["2"] = Q2;
  Q["3"] = Q3;
  Q["4"] = Q4;
  return(Q);
}


// [[Rcpp::export]]
List calcQ3_models(const arma::imat& M1,
	    const arma::imat& M2,
	    const arma::imat& M3, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
		   const NumericVector& pp3) { // pp for each model for disease 2
  // Rcpp::Rcout << "size of M: " << size(M1) << std::endl;
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  const int nmod3 = pp3.size();
  NumericMatrix Q1(nmod1,3);
  NumericMatrix Q2(nmod2,3);
  NumericMatrix Q3(nmod3,3);
  for(int i1=0; i1<nmod1; i1++) {
     for(int i2=0; i2<nmod2; i2++) {
       int idx1=modoverlap(M1.col(i1), M2.col(i2));
            for(int i3=0; i3<nmod3; i3++) {
   	int idx=idx1 +
   	  modoverlap(M1.col(i1), M3.col(i3)) +
   	  modoverlap(M2.col(i2), M3.col(i3));
	// if(idx > 2) {
   	// Rcpp::Rcout << i1 << ' ' << i2 << ' ' << i3 << ' ' << idx << std::endl;
	// }
	if(idx > 0) {
	  idx = idx - 1;
	  Q1(i1,idx) = Q1(i1,idx) + pp2(i2) * pp3(i3);
  	  Q2(i2,idx) = Q2(i2,idx) + pp1(i1) * pp3(i3);
  	  Q3(i3,idx) = Q3(i3,idx) + pp1(i1) * pp2(i2);
  	}
       }
     }
  }
  List Q;
  Q["1"] = Q1;
  Q["2"] = Q2;
  Q["3"] = Q3;
  return(Q);
}


// [[Rcpp::export]]
NumericVector calcQone2(const List S1, // model matrix - columns are models, rows are SNPs
	       const List S2, // model matrix - columns are models, rows are SNPs
	       const NumericVector& pp1,
	       const NumericVector& pp2) { // pp for each model for disease 2
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  // const int nsnp = M.nrow();
  NumericVector Q(nmod1);
  for(int i=0; i<nmod1; i++) {
    for(int j=0; j<nmod2; j++) {
      int idx=stroverlap( S1[i], S2[j]);
      if(idx > 0) {
	Q(i) = Q(i) + pp2(j);
      }
    }
  }
  return(Q);
}

// [[Rcpp::export]]
NumericMatrix calcQone3(const List S1,
	       const List S2,
	       const List S3, // model matrix - columns are models, rows are SNPs
	       const NumericVector& pp1,
	       const NumericVector& pp2,
	       const NumericVector& pp3) { // pp for each model for disease 2
  // Rcpp::Rcout << "size of M: " << size(M1) << std::endl;
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  const int nmod3 = pp3.size();
  NumericMatrix Q(nmod1,3);
  for(int i1=0; i1<nmod1; i1++) {
    for(int i2=0; i2<nmod2; i2++) {
      int idx1=stroverlap(S1[i1], S2[i2]) - 1;
      double pp12 = pp1(i1) * pp2(i2);
      for(int i3=0; i3<nmod3; i3++) {
   	int idx=idx1 +
   	  stroverlap(S1[i1], S3[i3]) +
   	  stroverlap(S2[i2], S3[i3]);
	if(idx > 0) 
	  Q(i1,idx) = Q(i1,idx) + pp2(i2) * pp3(i3);
      }
    }
  }
  return(Q);
}

// [[Rcpp::export]]
NumericMatrix calcQone4(const List S1,
	    const List S2,
	    const List S3, // model matrix - columns are models, rows are SNPs
	    const List S4, // model matrix - columns are models, rows are SNPs
	    const NumericVector& pp1,
	    const NumericVector& pp2,
	    const NumericVector& pp3,
	    const NumericVector& pp4) { // pp for each model for disease 2
  // Rcpp::Rcout << "size of M: " << size(M1) << std::endl;
  const int nmod1 = pp1.size();
  const int nmod2 = pp2.size();
  const int nmod3 = pp3.size();
  const int nmod4 = pp4.size();
  NumericMatrix Q(nmod1,6);
  for(int i1=0; i1<nmod1; i1++) {
     for(int i2=0; i2<nmod2; i2++) {
       int idx1=stroverlap(S1[i1], S2[i2]) - 1;
       double pp12 = pp1(i1) * pp2(i2);
       for(int i3=0; i3<nmod3; i3++) {
	 int idx2=idx1 +
	   stroverlap(S1[i1], S3[i3]) +
	   stroverlap(S2[i2], S3[i3]);
	 double pp23 = pp2(i2) * pp3(i3);
	 double pp13 = pp1(i1) * pp3(i3);
	 double pp123 = pp1(i1) * pp2(i2) * pp3(i3);
	 for(int i4=0; i4<nmod4; i4++) {
	   int idx=idx2 +
	     stroverlap(S1[i1], S4[i4]) +
	     stroverlap(S2[i2], S4[i4]) +
	     stroverlap(S3[i3], S4[i4]);
	   if(idx > 0) 
	     Q(i1,idx) = Q(i1,idx) + pp23 * pp4(i4);
	 }
       }
     }
  }
  return(Q);
}

*/
