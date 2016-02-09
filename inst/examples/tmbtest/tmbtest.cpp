
#include <TMB.hpp>
#include <covafill/TMB>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(coord);
  DATA_VECTOR(covObs);
  DATA_INTEGER(p);
  DATA_VECTOR(h);
  PARAMETER_VECTOR(x);
  covafill<Type> cf(coord,covObs,h,p);
  Type val = evalFill((CppAD::vector<Type>)x, cf)[0];
  return val;
}
