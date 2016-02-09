
#include <TMB.hpp>
#include <covafill/TMB>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(coord);
  DATA_VECTOR(covObs);
  DATA_INTEGER(p);
  DATA_VECTOR(h);
  DATA_SCALAR(d);
  PARAMETER_VECTOR(x);
  covafill<Type> cf(coord,covObs,h,p);
  covatree<Type> ct(d,&cf);
  Type val = evalTree((CppAD::vector<Type>)x, ct)[0];
  return val;
}
