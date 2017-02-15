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



#ifndef _COVAFILL_ATOMIC_TREE_
#define _COVAFILL_ATOMIC_TREE_

/** \brief Evaluates a covatree object, \a ct, at the coordinates \a tx.
 * \ingroup tmb
 */
CppAD::vector<double> evalTree(CppAD::vector<double> tx, const covatree<double> &ct)CSKIP({
    CppAD::vector<double> ty(1);
    // The first index from the operator is the function value
    ty[0] = ct.operator()(tx)[0];
    return ty;								
  })

/*! 
 * \overload
 * \ingroup tmb
 */
template <class Type>
CppAD::vector<Type> evalTree(CppAD::vector<Type> tx,const covatree<AD<Type> > &ct){
    CppAD::vector<Type> ty(1);
    // The first index from the operator is the function value
    ty[0] = CppAD::Value(ct.operator()(tx)[0]);
     return ty;								
  }



/** \brief CppAD atomic class to use estimated derivatives in automatic differentiation. See CppAD::atomic_base for further documentation.
 * \ingroup tmb
 */
template <class Type>							
class atomicEvalTree : public CppAD::atomic_base<Type> {		
public:
  /** \brief Constructs class to evaluate atomic function.*/
  atomicEvalTree(const char* name,covatree<AD<Type> > ct_) : CppAD::atomic_base<Type>(name), ct(ct_){
    atomic::atomicFunctionGenerated = true;				
    if(config.trace.atomic)						
    	std::cout << "Constructing atomic " << "evalTree" << "\n" ;	
    this->option(CppAD::atomic_base<Type>::bool_sparsity_enum);		
  }									
private:
  covatree<AD<Type> > ct;

  // Changed to use specific function value instead of using macro input
  virtual bool forward(size_t p,					
		       size_t q,					
		       const CppAD::vector<bool>& vx,			
		       CppAD::vector<bool>& vy,				
		       const CppAD::vector<Type>& tx,			
		       CppAD::vector<Type>& ty				
		       )						
  {									
    if(q>0)error("Atomic 'evalTree' order not implemented.\n");	
    if( vx.size() > 0 ){						
      bool anyvx = false;						
      for(size_t i=0;i<vx.size();i++)anyvx |= vx[i];			
      for(size_t i=0;i<vy.size();i++)vy[i] = anyvx;			
    }									
    ty = evalTree(tx, ct);	       					
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
    if(q>0)error("Atomic 'evalTree' order not implemented.\n");	
    CppAD::vector<AD<Type> > val = ct(tx);
    for(int i = 0; i < ct.getDim(); ++i)
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
CppAD::vector<AD<Type > > evalTree(CppAD::vector<AD<Type > > tx, covatree<AD<Type> > ct){
  static class atomicEvalTree<Type> afunEvalField("atomicEvalTree",ct); 
  CppAD::vector<AD<Type > > ty(1);				
  afunEvalField(tx,ty);						
  return ty;								
}


#endif
