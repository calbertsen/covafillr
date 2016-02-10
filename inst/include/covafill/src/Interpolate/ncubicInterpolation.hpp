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


template<typename scalartype_>
class ncubicInterpolation {
public:

  DEFINE_TYPES(scalartype_);

  ncubicInterpolation(covafill<scalartype>* cf,
		      vectortype minCoord_,
		      vectortype maxCoord_);

  virtual vectortype operator()(vectortype newcoord) = 0;

protected:
  ncubicInterpolation(vectortype minCoord_,
		      vectortype maxCoord_);

  int dim;
  vectortype minCoord;
  vectortype maxCoord;

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
};

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
};




#endif
