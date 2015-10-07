#ifndef _fill_ARR_H
#define _fill_ARR_H

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP fillArr(SEXP Ind_, SEXP yr_, SEXP dims_) {
  IntegerMatrix Ind(Ind_);
  NumericVector yr(yr_);
  IntegerVector dims(dims_);
  arma::cube myarr(dims[0],dims[1],dims[2]);
  for (unsigned int i = 0; i < yr.size(); i++) {
    myarr(Ind(i,0),Ind(i,1),Ind(i,2)) = yr[i];
  }
  return wrap(myarr);
}

#endif
