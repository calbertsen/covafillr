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

#ifndef _COVAFILL_NCUBIC_CLASS_
#define _COVAFILL_NCUBIC_CLASS_

/*! \brief Class for n-cubic interpolation (n = 1,2,3) of local polynomial regression on a square. The class should not be used as anything but a common parent for the dimension specific interpolation classes.
 *  \ingroup interpolate
*/
template<typename scalartype_>
class ncubicInterpolation {
  DEFINE_TYPES(scalartype_)


public:

  
  /** \brief Constructs a n-cubicInterpolation class from a covafill class \a cf, an boundaries of the interpolation square defined by the minimum coordinates, \a minCoord, and maximum coordinates, \a maxCoord, in each dimension, e.g., minCoord = (0,0) and maxCoord = (1,1). 
   */
  ncubicInterpolation(covafill<scalartype>* cf,
		      vectortype minCoord_,
		      vectortype maxCoord_);

  /** \brief Calculates the interpolation prediction at \a newcoord.
   */
  virtual vectortype operator()(vectortype newcoord) = 0;

protected:
  /** \brief Constructor from coordinates. Should in general not be called. */
  ncubicInterpolation(vectortype minCoord_,
		      vectortype maxCoord_);

  int dim; /**< Dimension of coordinates, i.e., the n in n-cubic. */
  vectortype minCoord; /**< Minimum coordinates of boundary box. */
  vectortype maxCoord; /**< maximum coordinates of boundary box. */

private:

  ncubicInterpolation<scalartype>* minChild;
  ncubicInterpolation<scalartype>* maxChild;



};


template<typename scalartype_>
ncubicInterpolation<scalartype_>::ncubicInterpolation(vectortype minCoord_,
						      vectortype maxCoord_)
  : dim(minCoord_.size()),
    minCoord(minCoord_),
    maxCoord(maxCoord_)
{
  minChild = NULL;
  maxChild = NULL;
}

  template<typename scalartype_>
ncubicInterpolation<scalartype_>::ncubicInterpolation(covafill<scalartype>* sf,
						      vectortype minCoord_,
						      vectortype maxCoord_)
  : dim(minCoord_.size()),
    minCoord(minCoord_),
    maxCoord(maxCoord_)
    
{
  minChild = NULL;
  maxChild = NULL;
}




#endif
