
#define CHECK_TREE_POINTER(name)		\
  if(R_ExternalPtrTag(name) != install("covatreePointer")) \
    Rf_error("The pointer must be to a covatree object"); \
  if(!R_ExternalPtrAddr(name)) \
    Rf_error("The pointer address is not valid");

extern "C" {

    static void finalizeTree(SEXP ptr){
    warning("Finalizing covatree");
    if(!R_ExternalPtrAddr(ptr))
      return;
    warning("Deleting covatree");
    delete (covatree<double>*)R_ExternalPtrAddr(ptr);
    R_ClearExternalPtr(ptr);    
  }
  
  SEXP MakeTree(SEXP coord,SEXP obs,SEXP h,SEXP p, SEXP d){
    covafill<double>* cfp = new covafill<double>(asMatrix(coord),
						       asVector(obs),
						       asVector(h),
						       asInteger(p));
    covatree<double>* ctp = new covatree<double>(asDouble(d),cfp);
    delete cfp;
    if(ctp == NULL){
      return R_NilValue;
    }
    SEXP val = R_MakeExternalPtr(ctp, install("covatreePointer"), R_NilValue);
    PROTECT(val);
    R_RegisterCFinalizerEx(val, finalizeTree, TRUE);
    UNPROTECT(1);
    return val;
  }

  SEXP getTreeDim(SEXP sp){
    CHECK_TREE_POINTER(sp);
    covatree<double>* ptr=(covatree<double>*)R_ExternalPtrAddr(sp);
    return(asSEXP(ptr->getDim()));
  }

  SEXP predictTree(SEXP sp, SEXP x){
    CHECK_TREE_POINTER(sp);
    covatree<double>* ptr=(covatree<double>*)R_ExternalPtrAddr(sp);
    int dim = ptr->getDim();
    
    if(isMatrix(x)){
      MatrixXd res(nrows(x),1 + dim);
      MatrixXd x0 = asMatrix(x);
      for(int i = 0; i < nrows(x); ++i)
	res.row(i) = ptr->operator()((vector)x0.row(i));
      return asSEXP(res);
    }else if(isNumeric(x)){
      return asSEXP(ptr->operator()(asVector(x)));
    }else{
      Rf_error("Element must be a matrix or numeric vector");
    }
    return R_NilValue;
  }
  
}
