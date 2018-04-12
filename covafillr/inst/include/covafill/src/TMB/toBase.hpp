// This file is part of covafill, a C++ template library for
// local polynomial regression of covariates in state-space models
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:

//     Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.

//     Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in
//     the documentation and/or other materials provided with the
//     distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef _COVAFILL_TMB_TO_BASE_
#define _COVAFILL_TMB_TO_BASE_


template<typename Type>
vector<Type> VectortoBase(vector<AD<Type> > x){
  vector<Type> res(x.size());
  for(int i = 0; i < res.size(); ++i)
    res(i) = CppAD::Value(x(i));
  return res;
}

vector<double> VectortoBase(vector<double> x){
  return x;
}

template<typename Type>
matrix<Type> MatrixtoBase(matrix<AD<Type> > x){
  matrix<Type> res(x.rows(), x.cols());
  for(int i = 0; i < res.rows(); ++i)
    for(int j = 0; j < res.cols(); ++j)
      res(i,j) = CppAD::Value(x(i,j));

  return res;
}

matrix<double> MatrixtoBase(matrix<double> x){
  return x;
}

template<typename Type>
covafill<Type> CFtoBase(const covafill<AD<Type> > &x){
  matrix<Type> coordinates = MatrixtoBase(x.coordinates);
  vector<Type> observations = VectortoBase(x.observations);
  int p = x.p;
  vector<Type> h = VectortoBase(x.h);
  return covafill<Type>(coordinates,observations,h,p);
}

covafill<double> CFtoBase(const covafill<double> &x){
  return x;
}

template<typename Type>
Type toBase(AD<Type> x){
  return CppAD::Value(x);
}

double toBase(double x){
  return x;
}



#endif
