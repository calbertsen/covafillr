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

#ifndef _COVAFILL_CUBIC_CLASS_
#define _COVAFILL_CUBIC_CLASS_

/** \defgroup interpolate Interpolation module
*
* The Interpolate module of covafill provides classes for cubic interpolation in 1-3 dimensions. The class serves as a common wrapper for the dimension specific interpolation classes.
* \verbatim
#include <covafill/Interpolate>
\endverbatim
*/

/*! \brief Class for cubic interpolation in dimension 1-3 of local polynomial regression on a square. 
 *  \ingroup interpolate
*/
template<typename scalartype_>
class cubicInterpolation {

  DEFINE_TYPES(scalartype_)


public:


  /** \brief Constructs a bicubicInterpolation class from a covafill class \a cf, an boundaries of the interpolation square defined by the minimum coordinates, \a minCoord, and maximum coordinates, \a maxCoord, in each dimension, e.g., minCoord = (0,0) and maxCoord = (1,1). 
   */
  cubicInterpolation(covafill<scalartype>* cf,
		     vectortype minCoords,
		     vectortype maxCoords);

  /** \brief Returns the interpolation prediction at \a newcoord.
   */
  
  vectortype operator()(vectortype newcoord);

private:


  ncubicInterpolation<scalartype>* nci; /**< Pointer to the n-cubic interpolation class which does all the work. */


};

template<typename scalartype_>
cubicInterpolation<scalartype_>::cubicInterpolation(covafill<scalartype>* cf,
						    vectortype minCoord,
						    vectortype maxCoord){
  int d = minCoord.size();
  if(d < 1){
    // Trow exception
    nci = NULL;
  }else if(d == 1){
    nci = new unicubicInterpolation<scalartype>(cf,minCoord,maxCoord);
  }else if(d == 2){
    nci = new bicubicInterpolation<scalartype>(cf,minCoord,maxCoord);
  }else if(d == 3){
    nci = new tricubicInterpolation<scalartype>(cf,minCoord,maxCoord);
  }else{
    nci = NULL;
  }

}

template<typename scalartype_>
typename cubicInterpolation<scalartype_>::vectortype cubicInterpolation<scalartype_>::operator()(vectortype newcoord) {
  return nci->operator()(newcoord);
}

#endif
