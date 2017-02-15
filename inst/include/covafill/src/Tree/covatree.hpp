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

#ifndef _COVATREE_BASE_CLASS_
#define _COVATREE_BASE_CLASS_

/** \defgroup tree Tree module
*
* The Tree module of covafill provides a class for search tree approximated local polynomial regression.
* \verbatim
#include <covafill/Interpolate>
\endverbatim
*/

/*! \brief Class that defines a covatree for search tree approximated local polynomial regression. 
 *  \ingroup tree
*/
template<typename scalartype_>
class covatree {

  DEFINE_TYPES(scalartype_)
  
public:

 
  // Constructors
  
  // Copy constructors
  // covatree(const covatree<scalartype_>& x);
  // #ifdef TMB_EXTERN
  // covatree(const covatree<AD<scalartype_> >& x);
  // #endif


  // Constructors that creates covafill first
  // covatree(matrixtype coordinates_,
  // 	     vectortype observations_,
  // 	     scalartype minSplitSize_);
  // covatree(matrixtype coordinates_,
  // 	     vectortype observations_,
  // 	     scalartype h_,
  // 	     int p_,
  // 	     scalartype minSplitSize_);
  // covatree(matrixtype coordinates_,
  // 	     vectortype observations_,
  // 	     vectortype h_,
  // 	     int p_,
  // 	     scalartype minSplitSize_);

  // Constructors from covafill pointer
  /** \brief Constructs a tree from a covafill object \a cf with minimum number of coordinates at which a sub tree will be created \a minSplitSize_. */
  covatree(scalartype minSplitSize_,
	     covafill<scalartype>* cf);

  /** \brief Destructor */
  ~covatree();
  
  // Public functions
  /** \brief Get coordinate dimension.  */
  int getDim();
  // getTree
  // getBoundingBoxes

  // Operators
  /** \brief Returns the interpolated value at \a newcoord of the local polynomial regressions at the corners of the boundary box.  */
  vectortype operator()(vectortype newcoord) const;
  // covatree<scalartype> & operator= (const covatree<scalartype>& rhs);

private:

  covanode<scalartype>* root;
  // lookup


};


template<typename scalartype_>
int covatree<scalartype_>::getDim() {
  return  root->getDim();
}


template<typename scalartype_>
typename covatree<scalartype_>::vectortype covatree<scalartype_>::operator()(vectortype newcoord) const {
  return  root->operator()(newcoord);
}





#endif
