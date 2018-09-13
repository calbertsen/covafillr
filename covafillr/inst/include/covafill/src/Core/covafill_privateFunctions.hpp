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

#ifndef _COVAFILL_PRIVATE_
#define _COVAFILL_PRIVATE_


template<typename scalartype_>
typename covafill<scalartype_>::scalartype covafill<scalartype_>::calcNorm(vectortype x0, vectortype x1) const {
  vectortype tmp = x0 - x1;
  scalartype ressq = 0.0;
  for(int i = 0; i < tmp.size(); ++i)
    ressq += tmp(i) * tmp(i);
  return sqrt(ressq);
}





template<typename scalartype_>
typename covafill<scalartype_>::scalartype covafill<scalartype_>::getWeight(vectortype x0, vectortype x1) const {

  scalartype_ sqnorm = (Hinv * (x0 - x1).matrix()).squaredNorm();

#ifdef USE_GAUSSIAN_KERNEL
  scalartype res = exp(-0.5 * sqnorm);
  #else
  scalartype a = 1.0 - sqnorm;
  scalartype b = 0.0;
  scalartype res = 0.5 * (a + b + fabs(a-b));
  #endif
  return res * detHinv;
}



#endif
