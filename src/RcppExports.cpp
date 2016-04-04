// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// vertboot_matrix_rcpp
IntegerMatrix vertboot_matrix_rcpp(IntegerMatrix m1, IntegerVector blist);
RcppExport SEXP snowboot_vertboot_matrix_rcpp(SEXP m1SEXP, SEXP blistSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type blist(blistSEXP);
    __result = Rcpp::wrap(vertboot_matrix_rcpp(m1, blist));
    return __result;
END_RCPP
}