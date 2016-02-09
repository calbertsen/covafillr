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

template<typename scalartype_>
class covafill {
  
public:
  DEFINE_TYPES(scalartype_);
  matrixtype coordinates;
  vectortype observations;
  int p;			// Degree
  vectortype h;

  
  // Constructors
  covafill(const covafill<scalartype_>& x);
  #ifdef TMB_EXTERN
  covafill(const covafill<AD<scalartype_> >& x);
  #endif
  covafill(matrixtype coordinates_,
	      vectortype observations_);
  covafill(matrixtype coordinates_,
	      vectortype observations_,
	      scalartype h_,
	      int p_);
  covafill(matrixtype coordinates_,
	      vectortype observations_,
	      vectortype h_,
	      int p_);

  // Public functions
  int getDim() const;
  void setH(scalartype h_);
  void setH(vectortype h_);


  // Operators
  vectortype operator()(vectortype x0, bool returnAll = false) const;
  vectortype operator()(vectortype x0, scalartype excludeRadius, bool returnAll = false) const;
  covafill<scalartype> & operator= (const covafill<scalartype>& rhs);

private:

  sparsematrixtype Hinv;	// Bandwiths
  scalartype detHinv;
  int dim;			// Dimension of coordinates (1, 2 or 3)
  int nobs;

  scalartype calcNorm(vectortype x0,
		      vectortype x1) const;
  scalartype getWeight(vectortype x0,
		       vectortype x1) const;

  
};

#endif
