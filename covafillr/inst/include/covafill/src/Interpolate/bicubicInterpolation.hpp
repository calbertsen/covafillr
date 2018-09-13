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

#ifndef _COVAFILL_BICUBIC_CLASS_
#define _COVAFILL_BICUBIC_CLASS_

/*! \brief Class for bi-cubic interpolation of local polynomial regression on a square. 
 *  \ingroup interpolate
*/
template<typename scalartype_>
class bicubicInterpolation : public ncubicInterpolation<scalartype_> {
  DEFINE_TYPES(scalartype_)

public:

  /** \brief Constructs a bicubicInterpolation class from a covafill class \a cf, an boundaries of the interpolation square defined by the minimum coordinates, \a minCoord, and maximum coordinates, \a maxCoord, in each dimension, e.g., minCoord = (0,0) and maxCoord = (1,1). 
   */
  bicubicInterpolation(covafill<scalartype>* cf,
		      vectortype minCoord,
		      vectortype maxCoord);

  virtual ~bicubicInterpolation();

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
typename bicubicInterpolation<scalartype_>::matrixtype bicubicInterpolation<scalartype_>::getMinv(){

  // Dim == 2
  //matrixtype Minv(16,16);
  scalartype Minv[16][16] =
    {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
     {-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
     {2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
     {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
     {0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
     {0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
     {-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
     {9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
     {-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
     {2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
     {-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
     {4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}};

  matrixtype MinvEigen(16,16);
  for(int i = 0; i < 16; ++i)
    for(int j = 0; j < 16; ++j)
      MinvEigen(i,j) = Minv[i][j];
			      
  return MinvEigen;
  
}

template<typename scalartype_>
typename bicubicInterpolation<scalartype_>::matrixtype bicubicInterpolation<scalartype_>::makeAlpha(covafill<scalartype>* cf,vectortype minCoord,vectortype maxCoord) {

  matrixtype coords(2,minCoord.size());
  for(int i = 0; i < minCoord.size(); ++i){
    coords(0,i) = minCoord(i);
    coords(1,i) = maxCoord(i);
  }
  
  matrixtype Minv = getMinv();

  int ncoef = 16;
  
  vectortype d = maxCoord - minCoord;
  vectortype mult(4);
  mult << 1, d(0), d(1), d(0) * d(1);
  vectortype x(ncoef);
  x.setZero();

  vectortype y(4), y1(4), y2(4), y12(4);
  int indxs[4][2] = {{0,0},{1,0},{1,1},{0,1}};
  
  for(int i = 0; i < 4; ++i){
    vectortype tmp(2);
    tmp(0) = coords(indxs[i][0],0);
    tmp(1) = coords(indxs[i][1],1);
    vectortype vals = cf->operator()(tmp, true);
    int nv = vals.size();
    vectortype tmp2(4);
    x(i) = vals(0)*mult(0); 	// f(x,y)
    x(i+4) = (nv > 1) ? vals(1)*mult(1) : 0.0;	// f_x(x,y)
    x(i+8) = (nv > 2) ? vals(2)*mult(2) : 0.0;	// f_y(x,y)
    x(i+12) = (nv > 5) ? vals(5)*mult(3) : 0.0;	// f_{x,y}(x,y)
  }
     

  vectortype alphaCoef = Minv * x.matrix();
  matrixtype coefs(4,ncoef/4);
  coefs.setZero();
  int indx = 0;
  for(int i = 0; i < 4; ++i)
    for(int j = 0; j < ncoef/4; ++j)
      coefs(i,j) = alphaCoef(indx++);

  return coefs;
}
  


template<typename scalartype_>
bicubicInterpolation<scalartype_>::bicubicInterpolation(covafill<scalartype>* cf,
						 vectortype minCoord,
						 vectortype maxCoord)
  : ncubicInterpolation<scalartype_>(minCoord,maxCoord),
    alpha(makeAlpha(cf,minCoord,maxCoord))
{}

template<typename scalartype_>
bicubicInterpolation<scalartype_>::~bicubicInterpolation(){
}

template<typename scalartype_>
typename bicubicInterpolation<scalartype_>::vectortype bicubicInterpolation<scalartype_>::operator()(vectortype newcoord) {
  vectortype newval(3);
  newval.setZero();

  vectortype d = this->maxCoord - this->minCoord;
  vectortype x = newcoord - this->minCoord;
  scalartype t = x(0)/d(0);
  scalartype u = x(1)/d(1);

  if(t > 1 || t < 0 || u > 1 || u < 0){
    // Trow exception
    return newval/(scalartype)0.0;
  }
    
  
  for(int i = 0; i < 4; ++i)
    for(int j = 0; j < 4; ++j){
       // Function value
      newval(0) += alpha(i,j) * pow(t,i) * pow(u,j);
      // Derivative of first variable - problem if t==0??
      if(i > 0)
	newval(1) += scalartype(i) * alpha(i,j) * pow(t,i-1) * pow(u,j) / d(0);
      // Derivative of second variable
      if(j > 0)
	newval(2) += scalartype(j) * alpha(i,j) * pow(t,i) * pow(u,j-1) / d(1);     
    }
  return newval;
}



#endif
