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

#ifndef _COVAFILL_OPERATORS_
#define _COVAFILL_OPERATORS_

#include<iostream>

template<typename scalartype_>
covafill<scalartype_> & covafill<scalartype_>::operator= (const covafill<scalartype_>& rhs){
  coordinates = rhs.coordinates;
  observations = rhs.observations;
  h = rhs.h;
  p = rhs.p;
  Hinv = rhs.Hinv;	// Bandwiths
  detHinv = rhs.detHinv;
  dim = rhs.dim;			// Dimension of coordinates (1, 2 or 3)
  nobs = rhs.nobs;
  return *this;
}



template<typename scalartype_>
typename covafill<scalartype_>::vectortype covafill<scalartype_>::operator()(vectortype x0, scalartype excludeRadius, bool returnAll) const {

  // Create (sparse) matrix of weights
  sparsematrixtype W(nobs,nobs);
  for(int i = 0; i < nobs; ++i){
    if(calcNorm(x0,coordinates.row(i)) > excludeRadius){
      scalartype wi = getWeight(x0,coordinates.row(i));
      if(wi > 0)
	W.insert(i,i) = sqrt(wi);	// = wi
    }
  }

  // Calculate (dense) design matrix for least squares (nobs x dim)
  int lsdim = 1 + dim;
  if(p >= 2)
    lsdim += 0.5 * dim * (dim + 1);
  if(p >= 3)
    lsdim += (p - 2) * dim;
  matrixtype X(nobs, lsdim);

  for(int i = 0; i < nobs; ++i){
    vectortype xtmp = (vectortype)coordinates.row(i) - x0;
    // First is the intercept
    X(i,0) = 1.0;
    // Then (x-x0)
    for(int j = 1; j < 1+dim; ++j)
      X(i,j) = xtmp(j-1);
    // Finally (x-x0)^T(x-x0) if p == 2; i.e. the interactions
    if(p >= 2){
      int indx = 1+dim;
      for(int j = 0; j < dim; ++j)
	for(int k = j; k < dim; ++k){
	  X(i,indx++) = xtmp(j) * xtmp(k) / 2.0;
	}
      if(p > 2){
	int fac = 2;
	for(int k = 3; k <= p; ++k){
	  fac *= k;
	  for(int j = 0; j < dim; ++j)
	    X(i,indx++) = pow(xtmp(j),k) / fac;
	}
      }
    }
  }

  // Calculate weighted least squares estimates
  matrixtype XTW = X.transpose() * W;
  matrixtype XWX = XTW * X;
  matrixtype tmp = XWX.inverse() * XTW;
  vectortype beta = tmp * observations.matrix();
  // Return the vector (f(x0), \partial/\partial x_1 f(x0), ..., \partial/\partial x_dim)

  if(returnAll)
    return beta;
  return beta.segment(0,1+dim);
};

template<typename scalartype_>
typename covafill<scalartype_>::vectortype covafill<scalartype_>::operator()(vectortype x0, bool returnAll) const {
  return operator()(x0,scalartype(-1.0),returnAll);
}

#endif