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

// This file is inspired by the definition of the TMB macro: TMB_ATOMIC_VECTOR_FUNCTION
// https://github.com/kaskr/adcomp/blob/f152ccae08b6720d86e43fe659a4a4067322368d/TMB/inst/include/atomic_macro.hpp


#ifndef _COVAFILL_ATOMIC_
#define _COVAFILL_ATOMIC_


/** \defgroup tmb TMB module
*
* The TMB module of covafill provides functions to evaluate a covafill or covatree object from a TMB model such that the estimated gradients are used in the automatic differentiation.
* \verbatim
#include <covafill/TMB>
\endverbatim
*/


/** \brief Evaluates a covafill object, \a cf, at the coordinates \a tx.
 * \ingroup tmb
 */
CppAD::vector<double> evalFill(CppAD::vector<double> tx, const covafill<double> &cf)CSKIP({
    CppAD::vector<double> ty(1);
    // The first index from the operator is the function value
    ty[0] = cf.operator()(tx)[0];
    return ty;								
  })

/*! 
 * \overload
 * \ingroup tmb
 */
template <class Type>
CppAD::vector<Type> evalFill(CppAD::vector<Type> tx,const covafill<AD<Type> > &cf){
    CppAD::vector<Type> ty(1);
    // The first index from the operator is the function value
    ty[0] = CppAD::Value(cf.operator()(tx)[0]);
     return ty;								
}
 

/** \brief CppAD atomic class to use estimated derivatives in automatic differentiation. See CppAD::atomic_base for further documentation.
 * \ingroup tmb
 */
template <class Type>							
class atomicEvalFill : public CppAD::atomic_base<Type> {		
public:

  /** \brief Constructs class to evaluate atomic function.*/
  atomicEvalFill(const char* name,covafill<AD<Type> > cf_) : CppAD::atomic_base<Type>(name), cf(cf_){
    atomic::atomicFunctionGenerated = true;				
    if(config.trace.atomic)						
      std::cout << "Constructing atomic " << "evalFill" << "\n" ;	
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);		
  }
 
private:
  covafill<AD<Type> > cf;

  // Changed to use specific function value instead of using macro input
  virtual bool forward(size_t p,					
		       size_t q,					
		       const CppAD::vector<bool>& vx,			
		       CppAD::vector<bool>& vy,				
		       const CppAD::vector<Type>& tx,			
		       CppAD::vector<Type>& ty				
		       )						
  {									
    if(q>0)error("Atomic 'evalFill' order not implemented.\n");	
    if( vx.size() > 0 ){						
      bool anyvx = false;						
      for(size_t i=0;i<vx.size();i++)anyvx |= vx[i];			
      for(size_t i=0;i<vy.size();i++)vy[i] = anyvx;			
    }									
    ty = evalFill(tx, cf);	       					
    return true;							
  }

  // Changed to calculate specific gradient instead of using macro input
  virtual bool reverse(size_t q,					
		       const CppAD::vector<Type>& tx,			
		       const CppAD::vector<Type>& ty,			
		       CppAD::vector<Type>& px,				
		       const CppAD::vector<Type>& py			
		       )						
  {									
    if(q>0)error("Atomic 'evalFill' order not implemented.\n");	
    CppAD::vector<AD<Type> > val = cf(tx);
    for(int i = 0; i < cf.getDim(); ++i)
      px[i] = CppAD::Value(val[i+1]) * py[0];
    return true;							
  }

  // Not changed from TMB
  virtual bool rev_sparse_jac(size_t q,					
			      const CppAD::vector<bool>& rt,		
			      CppAD::vector<bool>& st)			
  {									
    bool anyrt = false;							
    for(size_t i=0;i<rt.size();i++)anyrt |= rt[i];			
    for(size_t i=0;i<st.size();i++)st[i]=anyrt;				
    return true;							
  }

  // Not changed from TMB
  virtual bool rev_sparse_jac(size_t q,					
			      const CppAD::vector< std::set<size_t> >& rt, 
			      CppAD::vector< std::set<size_t> >& st)	
  {									
    error("Should not be called");					
  }									
};

/*! 
 * \overload
 * \ingroup tmb
 */
template<class Type> 							
CppAD::vector<AD<Type > > evalFill(CppAD::vector<AD<Type > > tx,covafill<AD<Type> > cf){
  static class atomicEvalFill<Type> afunEvalFill("atomicEvalFill",cf); 
  CppAD::vector<AD<Type > > ty(1);				
  afunEvalFill(tx,ty);						
  return ty;								
}


#endif
