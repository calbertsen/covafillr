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

#ifndef _COVAFILL_BASE_CLASS_
#define _COVAFILL_BASE_CLASS_

/** \defgroup core Core module
*
* The Core module of covafill provides a class for local polynomial regression.
* \verbatim
#include <covafill/Core>
\endverbatim
*/

/*! \brief Class to do local polynomial regression. 
 *  \ingroup core
*/
template<typename scalartype_>
class covafill {
  DEFINE_TYPES(scalartype_)
public:
  matrixtype coordinates;  /**< Coordinates/covariates of input. */
  vectortype observations; /**< Input observations. */
  int p;                   /**< Polynomial degree. */
  vectortype h;            /**< Vector of (positive) bandwiths - one for each covariate.*/

  
  // Constructors

  /** \brief Constructs a covafill class from another covafill class \a x. 
   */
  covafill(const covafill<scalartype_>& x);
  #ifdef TMB_EXTERN
  covafill(const covafill<AD<scalartype_> >& x);
  #endif

  /** \brief Constructs a covafill class with coordinates matrix \a coordinates_, observation vector \a obervations, bandwiths 1, and polynomial degree 2. 
   */
  covafill(matrixtype coordinates_,
	      vectortype observations_);

  /** \brief Constructs a covafill class with coordinates matrix \a coordinates_, observation vector \a obervations, bandwiths \a h_, and polynomial degree \a p_. 
   */
  covafill(matrixtype coordinates_,
	      vectortype observations_,
	      scalartype h_,
	      int p_);

  /** \brief Constructs a covafill class with coordinates matrix \a coordinates_, observation vector \a obervations, bandwiths \a h_, and polynomial degree \a p_. 
   */
  covafill(matrixtype coordinates_,
	      vectortype observations_,
	      vectortype h_,
	      int p_);

  // Public functions
  int getDim() const;          /*!< Returns the covariate dimension. */
  void setH(scalartype h_);    /*!< Sets all bandwiths to h_. */
  void setH(vectortype h_);    /*!< Sets the bandwiths from a vector. The length of h_ must match the covariate dimension. */


  // Operators

  /** \brief Calculates the local polynomial regression estimate at \a x0. If \a returnAll is false, then only the function and first derivative estimates are returned. Otherwise all estimates are returned. */
  vectortype operator()(vectortype x0, bool returnAll = false) const;

  /** \brief Calculates the local polynomial regression estimate and covariance of estimates at \a x0. If \a returnAll is false, then only the function and first derivative estimates are returned. Otherwise all estimates are returned. For cov = 0, only the variance of the estimates are returned, for cov = 1 the covariance matrix is returned. If returnAll is false, cov is ignored and asumed to be 0.*/
  vecvectype operator()(vectortype x0, int cov, bool returnAll = false) const;
  
  /** \brief Calculates the local polynomial regression estimate at \a x0. All observations with coordinates \f$ x \f$ such that \f$ \|x-x_0\| > r \f$, where \a r is \a exludeRadius. If \a returnAll is false, then only the function and first derivative estimates are returned. Otherwise all estimates are returned. */
  vectortype operator()(vectortype x0, scalartype excludeRadius, bool returnAll = false) const;

  /** \brief Assignment operator for covafill.*/
  covafill<scalartype> & operator= (const covafill<scalartype>& rhs);

private:
  sparsematrixtype Hinv;	/**< Matrix of inverse bandwiths. */
  scalartype detHinv;           /**< Determinant of Hinv. */
  int dim;			/**< Number of columns of \a coordinates. */
  int nobs;                     /**< Number of observations, i.e., length og \a observations. */

  /** \brief Calculates \f$ \|x_0-x_1\|_2 \f$. */
  scalartype calcNorm(vectortype x0,
		      vectortype x1) const;
  /** \brief Calculates the multivariate epanechnikov kernel between \a x0 and \a x1.
   */
  scalartype getWeight(vectortype x0,
		       vectortype x1) const;

  
};

#endif
