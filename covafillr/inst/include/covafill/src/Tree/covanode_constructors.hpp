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

#ifndef _COVANODE_CONSTRUCTORS_
#define _COVANODE_CONSTRUCTORS_

template<typename scalartype_>
covanode<scalartype_>::covanode(matrixtype coordSplit,
				    scalartype minSplitSize_,
				    covafill<scalartype>* cf,
				    vectortype minCoords,
				    vectortype maxCoords) : minSplitSize(minSplitSize_), dim(minCoords.size()) {

  if(coordSplit.rows() < minSplitSize_){
     // TERMINAL NODE
    isTerminal = bool(true);
    left = NULL;
    right = NULL;
    cubic = new cubicInterpolation<scalartype>(cf,minCoords,maxCoords);
  }else{
    // ROOT OF SUBTREE
    isTerminal = bool(false);
    vectortype xx(coordSplit.cols());
    xx.setZero();
    vectortype xx2(coordSplit.cols());
    xx2.setZero();
    for(int i = 0; i < xx.size(); ++i)
      for(int j = 0; j < coordSplit.rows(); ++j){
	xx(i) += coordSplit(j,i);
	xx2(i) += coordSplit(j,i) * coordSplit(j,i);
      }
    vectortype mn = xx/(scalartype)coordSplit.rows(); // Change to median?
    vectortype vr = xx2/(scalartype)coordSplit.rows() - mn * mn;
    int d_ = 0;
    for(int i = 1; i < vr.size(); ++i)
      if(vr(i) > vr(d_))
	d_ = i;

    splitDim = int(d_);
    splitVal = scalartype(mn(d_));
 
    // Split data in left and right
    int lngth[2] = {0,0};
    
    for(int i = 0; i < coordSplit.rows(); ++i)
      (coordSplit(i,splitDim) <= splitVal) ? lngth[0]++ : lngth[1]++;

    matrixtype splitLeft(lngth[0],coordSplit.cols());
    matrixtype splitRight(lngth[1],coordSplit.cols());

    int lindx = 0;
    int rindx = 0;
    
    for(int i = 0; i < coordSplit.rows(); ++i)
      if(coordSplit(i,splitDim) <= splitVal){
	splitLeft.row(lindx++) = coordSplit.row(i);
      }else{
	splitRight.row(rindx++) = coordSplit.row(i);
      }
    vectortype lmin = minCoords;
    vectortype rmin = minCoords;
    vectortype lmax = maxCoords;
    vectortype rmax = maxCoords;

    // For the left split, the d'th variable is now less than splitVal
    lmax(splitDim) = splitVal;
    // For the right split, the d'th variable is now greater than splitVal
    rmin(splitDim) = splitVal;

    // Recursion
    left = new covanode(splitLeft,minSplitSize_,cf,lmin,lmax);
    right = new covanode(splitRight,minSplitSize_,cf,rmin,rmax);
    cubic = NULL;
  }

}


template<typename scalartype_>
covanode<scalartype_>::~covanode(){
  delete left;
  delete right;
  delete cubic;
}

#endif
