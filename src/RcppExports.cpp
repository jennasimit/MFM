// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// modoverlap
int modoverlap(const arma::ivec& x, const arma::ivec& y);
RcppExport SEXP _MTFM_modoverlap(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(modoverlap(x, y));
    return rcpp_result_gen;
END_RCPP
}
// stroverlap
int stroverlap(const IntegerVector& x, const IntegerVector& y);
RcppExport SEXP _MTFM_stroverlap(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(stroverlap(x, y));
    return rcpp_result_gen;
END_RCPP
}
// strint
int strint(const IntegerVector& x, const IntegerVector& y);
RcppExport SEXP _MTFM_strint(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(strint(x, y));
    return rcpp_result_gen;
END_RCPP
}
// calcQ2
List calcQ2(const List S1, const List S2, const NumericVector& pp1, const NumericVector& pp2);
RcppExport SEXP _MTFM_calcQ2(SEXP S1SEXP, SEXP S2SEXP, SEXP pp1SEXP, SEXP pp2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< const List >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp1(pp1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp2(pp2SEXP);
    rcpp_result_gen = Rcpp::wrap(calcQ2(S1, S2, pp1, pp2));
    return rcpp_result_gen;
END_RCPP
}
// calcQ2_models
List calcQ2_models(const arma::imat& M1, const arma::imat& M2, const NumericVector& pp1, const NumericVector& pp2);
RcppExport SEXP _MTFM_calcQ2_models(SEXP M1SEXP, SEXP M2SEXP, SEXP pp1SEXP, SEXP pp2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type M2(M2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp1(pp1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp2(pp2SEXP);
    rcpp_result_gen = Rcpp::wrap(calcQ2_models(M1, M2, pp1, pp2));
    return rcpp_result_gen;
END_RCPP
}
// calcQ3
List calcQ3(const List S1, const List S2, const List S3, const NumericVector& pp1, const NumericVector& pp2, const NumericVector& pp3);
RcppExport SEXP _MTFM_calcQ3(SEXP S1SEXP, SEXP S2SEXP, SEXP S3SEXP, SEXP pp1SEXP, SEXP pp2SEXP, SEXP pp3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< const List >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const List >::type S3(S3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp1(pp1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp2(pp2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp3(pp3SEXP);
    rcpp_result_gen = Rcpp::wrap(calcQ3(S1, S2, S3, pp1, pp2, pp3));
    return rcpp_result_gen;
END_RCPP
}
// calcQ3log
List calcQ3log(const List S1, const List S2, const List S3, const NumericVector& pp1, const NumericVector& pp2, const NumericVector& pp3);
RcppExport SEXP _MTFM_calcQ3log(SEXP S1SEXP, SEXP S2SEXP, SEXP S3SEXP, SEXP pp1SEXP, SEXP pp2SEXP, SEXP pp3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< const List >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const List >::type S3(S3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp1(pp1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp2(pp2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp3(pp3SEXP);
    rcpp_result_gen = Rcpp::wrap(calcQ3log(S1, S2, S3, pp1, pp2, pp3));
    return rcpp_result_gen;
END_RCPP
}
// calcQ4
List calcQ4(const List S1, const List S2, const List S3, const List S4, const NumericVector& pp1, const NumericVector& pp2, const NumericVector& pp3, const NumericVector& pp4);
RcppExport SEXP _MTFM_calcQ4(SEXP S1SEXP, SEXP S2SEXP, SEXP S3SEXP, SEXP S4SEXP, SEXP pp1SEXP, SEXP pp2SEXP, SEXP pp3SEXP, SEXP pp4SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< const List >::type S2(S2SEXP);
    Rcpp::traits::input_parameter< const List >::type S3(S3SEXP);
    Rcpp::traits::input_parameter< const List >::type S4(S4SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp1(pp1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp2(pp2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp3(pp3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp4(pp4SEXP);
    rcpp_result_gen = Rcpp::wrap(calcQ4(S1, S2, S3, S4, pp1, pp2, pp3, pp4));
    return rcpp_result_gen;
END_RCPP
}
// calcQ3_models
List calcQ3_models(const arma::imat& M1, const arma::imat& M2, const arma::imat& M3, const NumericVector& pp1, const NumericVector& pp2, const NumericVector& pp3);
RcppExport SEXP _MTFM_calcQ3_models(SEXP M1SEXP, SEXP M2SEXP, SEXP M3SEXP, SEXP pp1SEXP, SEXP pp2SEXP, SEXP pp3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type M1(M1SEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type M2(M2SEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type M3(M3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp1(pp1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp2(pp2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pp3(pp3SEXP);
    rcpp_result_gen = Rcpp::wrap(calcQ3_models(M1, M2, M3, pp1, pp2, pp3));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MTFM_modoverlap", (DL_FUNC) &_MTFM_modoverlap, 2},
    {"_MTFM_stroverlap", (DL_FUNC) &_MTFM_stroverlap, 2},
    {"_MTFM_strint", (DL_FUNC) &_MTFM_strint, 2},
    {"_MTFM_calcQ2", (DL_FUNC) &_MTFM_calcQ2, 4},
    {"_MTFM_calcQ2_models", (DL_FUNC) &_MTFM_calcQ2_models, 4},
    {"_MTFM_calcQ3", (DL_FUNC) &_MTFM_calcQ3, 6},
    {"_MTFM_calcQ3log", (DL_FUNC) &_MTFM_calcQ3log, 6},
    {"_MTFM_calcQ4", (DL_FUNC) &_MTFM_calcQ4, 8},
    {"_MTFM_calcQ3_models", (DL_FUNC) &_MTFM_calcQ3_models, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_MTFM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
