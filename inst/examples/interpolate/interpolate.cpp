
#include <TMB.hpp>
#include <covafill/TMB>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(coord);
  DATA_INTEGER(p);
  DATA_VECTOR(h);
  DATA_SCALAR(d);
  DATA_VECTOR(obs);
  DATA_MATRIX(x);
  DATA_VECTOR(repx);
  DATA_VECTOR(repy);  
  PARAMETER_VECTOR(covObs);
  PARAMETER(logSd);
  covafill<Type> cf(coord,covObs,h,p);
  covatree<Type> ct(d,&cf);
  Type nll = 0.0;
  for(int i = 0; i < obs.size(); ++i){
    vector<Type> xtmp = x.col(i);
    vector<Type> ctOut = ct(xtmp);
    nll -= dnorm(obs(i),ctOut(0),exp(logSd),true);
  }

  matrix<Type> predXY(repx.size(),repy.size());
  predXY.setZero();
  for(int i = 0; i < repx.size(); ++i)
    for(int j = 0; j < repx.size(); ++j){
      vector<Type> tmp(2);
      tmp(0) = repx(i); tmp(1) = repy(j);
      predXY(i,j) = cf(tmp)(0);
    }

  REPORT(predXY);
  return nll;
}
