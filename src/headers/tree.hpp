
extern "C" {
  SEXP MakeTree(SEXP coord,SEXP obs,SEXP h,SEXP p, SEXP d){


    covafill<double>* cfp = new covafill<double>(asMatrix(coord),
						       asVector(obs),
						       asVector(h),
						       asInteger(p));
    covatree<double>* ctp = new covatree<double>(asDouble(d),cfp);

    if(ctp == NULL){
      return R_NilValue;
    }
    SEXP val = R_MakeExternalPtr(ctp, install("covatreePointer"), R_NilValue);
    return val;
  }

  SEXP getTreeDim(SEXP sp){
    if(R_ExternalPtrTag(sp) != install("covatreePointer"))
      Rf_error("The pointer must be to a covatree object");
    covatree<double>* ptr=(covatree<double>*)R_ExternalPtrAddr(sp);
    return(asSEXP(ptr->getDim()));
  }

  SEXP predictTree(SEXP sp, SEXP x){
    if(R_ExternalPtrTag(sp) != install("covatreePointer"))
      Rf_error("The pointer must be to a covatree object");
 
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
