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

#ifndef _COVAFILL_UNICUBIC_CLASS_
#define _COVAFILL_UNICUBIC_CLASS_

/*! \brief Class for cubic interpolation of local polynomial regression on a square. 
 *  \ingroup interpolate
*/
template<typename scalartype_>
class unicubicInterpolation : public ncubicInterpolation<scalartype_> {

  DEFINE_TYPES(scalartype_)

public:


  /** \brief Constructs a unicubicInterpolation class from a covafill class \a cf, an boundaries of the interpolation square defined by the minimum coordinates, \a minCoord, and maximum coordinates, \a maxCoord, in each dimension, e.g., minCoord = 0 and maxCoord = 1. 
   */
  unicubicInterpolation(covafill<scalartype>* cf,
		      vectortype minCoord,
		      vectortype maxCoord);

  virtual ~unicubicInterpolation();
  /** \brief Calculates the interpolation prediction at \a newcoord.
   */
  virtual vectortype operator()(vectortype newcoord);

  
private:

  /*
    From ncubicInterpolation:
    int dim;
    vectortype minCoord;
    vectortype maxCoord;
  */
  matrixtype alpha; /**< Matrix of interpolation coordinates. */

  /** Returns \f$ M^{-1} \f$ (See Lekien & Marsden  2005) */
  matrixtype getMinv();
  
  /** Calculates interpolation coefficients (See Lekien & Marsden  2005) */
  matrixtype makeAlpha(covafill<scalartype>* cf,
		      vectortype minCoord,
		      vectortype maxCoord);

  
};


template<typename scalartype_>
typename unicubicInterpolation<scalartype_>::matrixtype unicubicInterpolation<scalartype_>::getMinv(){

  // Dim == 1
  // matrixtype Minv(4,4);
  // Minv << 1, 0, 0, 0,
  //   0, 1, 0, 0,
  //   -3, -2, 3, -1,
  //   2, 1, -2, 1;

  scalartype Minv[4][4] =
    {{ 1, 0, 0, 0},
     {0, 1, 0, 0},
     {-3, -2, 3, -1},
     {2, 1, -2, 1}};
  
  matrixtype MinvEigen(4,4);
  for(int i = 0; i < 4; ++i)
    for(int j = 0; j < 4; ++j)
      MinvEigen(i,j) = Minv[i][j];

  return MinvEigen;
}

template<typename scalartype_>
typename unicubicInterpolation<scalartype_>::matrixtype unicubicInterpolation<scalartype_>::makeAlpha(covafill<scalartype>* cf,vectortype minCoord,vectortype maxCoord) {
 
  matrixtype coords(2,minCoord.size());
  for(int i = 0; i < minCoord.size(); ++i){
    coords(0,i) = minCoord(i);
    coords(1,i) = maxCoord(i);
  }
  
  matrixtype Minv = getMinv();

  int ncoef = 4;
  
  vectortype d = maxCoord - minCoord;
  vectortype mult(2);
  mult << 1, d(0);

  vectortype x(ncoef);
  x.setZero();
  
  for(int i = 0; i < coords.rows(); ++i){
    vectortype tmp(1);
    tmp(0) = coords(i,0);
    vectortype vals = cf->operator()(tmp);
    x.segment(2*i,2) = vals.head(2) * mult;
  }

  vectortype alphaCoef = Minv * x.matrix();
  matrixtype coefs(4,ncoef/4);
  coefs.setZero();
  int indx = 0;
  for(int j = 0; j < ncoef/4; ++j)
    for(int i = 0; i < 4; ++i)
      coefs(i,j) = alphaCoef(indx++);

  return Minv;
}
  


template<typename scalartype_>
unicubicInterpolation<scalartype_>::unicubicInterpolation(covafill<scalartype>* cf,
						 vectortype minCoord,
						 vectortype maxCoord)
  : ncubicInterpolation<scalartype_>(minCoord,maxCoord),
    alpha(makeAlpha(cf,minCoord,maxCoord))
{}

template<typename scalartype_>
unicubicInterpolation<scalartype_>::~unicubicInterpolation(){
}


template<typename scalartype_>
typename unicubicInterpolation<scalartype_>::vectortype unicubicInterpolation<scalartype_>::operator()(vectortype newcoord) {
  
  vectortype newval(2);
  newval.setZero();

  vectortype d = this->maxCoord - this->minCoord;
  vectortype x = newcoord - this->minCoord;
  scalartype t = x(0)/d(0);

  if(t > 1 || t < 0){
    // Trow exception
    return newval/(scalartype)0.0;
  }

  
  for(int i = 0; i < 4; ++i){
    // Function value
    newval(0) += alpha(i,0) * pow(t,i);
    // Derivative of first variable - problem if t==0??
    if(i > 0)
      newval(1) += i * alpha(i,0) * pow(t,i-1) / d(0);
    // Derivative of second variable
  }
  return newval;
}


#endif
