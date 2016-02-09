
#ifndef _COVAFILL_R_UTILS_CONVERT_
#define _COVAFILL_R_UTILS_CONVERT_



#include <Eigen/Dense>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


// Inspired by or copied from github.com/kaskr/adcomp


using namespace Eigen;
typedef Array<double,Dynamic,1> vector;

SEXP asSEXP(double x){
  SEXP y = PROTECT(allocVector(REALSXP,1));
  REAL(y)[0] = x;

  UNPROTECT(1);
  return y; 

}

SEXP asSEXP(int x){
  SEXP y = PROTECT(allocVector(INTSXP,1));
  INTEGER(y)[0] = x;

  UNPROTECT(1);
  return y; 
}

SEXP asSEXP(vector x){
  int n = x.size();
  SEXP y = PROTECT(allocVector(REALSXP,n));

  for(int i = 0; i < n; ++i)
    REAL(y)[i] = x[i];

  UNPROTECT(1);
  return y; 
}

SEXP asSEXP(MatrixXd x){
  int nr = x.rows();
  int nc = x.cols();
  SEXP y = PROTECT(allocMatrix(REALSXP,nr,nc));
  for(int i = 0; i < nr; ++i)
    for(int j = 0; j < nc; ++j)
      REAL(y)[i+nr*j] = x(i,j);

  UNPROTECT(1);
  return y;
}

double asDouble(SEXP x){
  if(!(isNumeric(x) && length(x)==1))
    Rf_error("Element must be a numeric of length 1");
  return REAL(x)[0];
}

int asInteger(SEXP x){
  if(!(isInteger(x) && length(x)==1))
    Rf_error("Element must be an integer vector of length 1");
  return INTEGER(x)[0];
}

vector asVector(SEXP x) {

  if(!isNumeric(x))
    Rf_error("Element must be a numeric vector");

  int n = length(x);
  vector y(n);

  for(int i = 0; i < n; ++i)
    y(i) = REAL(x)[i];

  return y;
}

MatrixXd asMatrix(SEXP x)
{

  if(!isMatrix(x))
    Rf_error("Element must be a matrix");

  int nr = nrows(x);
  int nc = ncols(x);
  MatrixXd y(nr,nc);

  for(int i = 0; i < nr; ++i)
    for(int j = 0; j < nc; ++j)
      y(i,j)=REAL(x)[i+nr*j];

  return y;
}

#endif
