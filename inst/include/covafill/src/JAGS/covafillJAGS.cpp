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

#ifndef _COVAFILL_JAGS_FUN_
#define _COVAFILL_JAGS_FUN_

using std::vector;

typedef Eigen::Array<double,Eigen::Dynamic,1>  cVector;
typedef Eigen::MatrixXd cMatrix;

namespace jags {
  namespace covafillJAGS {

    covafillJAGS::covafillJAGS()
      : ArrayFunction("covafill", 5)
    {
    }

    void 
    covafillJAGS::evaluate (double *value, vector<double const *> const &args,
		       vector<vector<unsigned int> > const &dims) const
    {
   
      // Create X
      unsigned int nrowX = dims[0][0];
      unsigned int ncolX = dims[0].size() == 2 ? dims[0][1] : 1;
      cMatrix X(nrowX,ncolX);
      unsigned int indx;
      indx = 0;
      // JAGS uses column-major
      for (unsigned int j = 0; j < ncolX; ++j) {
	for (unsigned int i = 0; i < nrowX; ++i) {
	  X(i,j) = args[0][indx++];
	}
      }
      
	// Create coordinates
	unsigned int nrowC = dims[1][0];
	unsigned int ncolC = dims[1].size() == 2 ? dims[1][1] : 1;
	cMatrix coord(nrowC,ncolC);
	indx = 0;
	// JAGS uses column-major
	for (unsigned int j = 0; j < ncolC; ++j) {
	  for (unsigned int i = 0; i < nrowC; ++i) {
	    coord(i,j) = args[1][indx++];
	  }
	}
  

      // Create obs
      unsigned int lengthO = dims[2][0];
      cVector obs(lengthO);
      for(unsigned int i = 0; i < lengthO; ++i){
	obs[i] = args[2][i];
      }

      // Create h
      unsigned int lengthH = dims[3][0];
      cVector h(ncolC);
      for(unsigned int i = 0; i < ncolC; ++i){
	if(lengthH == 1){
	  h[i] = args[3][0];
	}else{
	  h[i] = args[3][i];
	}
      }

      // Create p
      double p = args[4][0];

      // Create covarfill
      covafill<double> cf(coord,obs,h,p);

      // Calculate
      for(unsigned int i = 0; i < nrowX; ++i)
	value[i] = cf((cVector)X.row(i))(0);      
    }

    vector<unsigned int> 
    covafillJAGS::dim (vector <vector<unsigned int> > const &dims,
		  vector<double const *> const &values) const
    {
      vector<unsigned int> ans(1);
      ans[0] = dims[0][0];

      return drop(ans);
    }

    bool 
    covafillJAGS::checkParameterDim (vector<vector<unsigned int> > const &dims) const
    {
      return true;// (isMatrix(dims[0]) || isVector(dims[0]))
	// && (isMatrix(dims[1]) || isVector(dims[1]))
	// && isVector(dims[2])
	// && (isVector(dims[3]) || isScalar(dims[3]))
	// && isScalar(dims[4]);
    } 
  }
}

#endif
