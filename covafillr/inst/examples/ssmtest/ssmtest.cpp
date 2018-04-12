
#include <TMB.hpp>
#include <covafill/TMB>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(obs);
  DATA_MATRIX(coord);
  DATA_VECTOR(covObs);
  DATA_INTEGER(p);
  DATA_VECTOR(h);

  PARAMETER(logObsSd);
  PARAMETER(logObsTSd);
  PARAMETER(logStatSd);
  PARAMETER_MATRIX(x);

  Type nll = 0.0;
  covafill<Type> cf(coord,covObs,h,p);

  // Contribution from states
  for(int i = 1; i < x.cols(); ++i){
    nll -= dnorm(x(0,i), x(0,i-1), exp(logStatSd),true);
    nll -= dnorm(x(1,i), x(1,i-1), exp(logStatSd),true);
  }

  // contribution from observations
  for(int i = 0; i < obs.cols(); ++i){
    nll -= dnorm(obs(0,i), x(0,i), exp(logObsSd),true);
    nll -= dnorm(obs(1,i), x(1,i), exp(logObsSd),true);
    vector<Type> tmp = x.col(i);
    Type val = evalFill((CppAD::vector<Type>)tmp, cf)[0];
    nll -= dnorm(obs(2,i), val, exp(logObsTSd),true);
  }
 
    
  return nll;
}
