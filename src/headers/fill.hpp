
extern "C" {
  SEXP MakeFill(SEXP coord,SEXP obs,SEXP h,SEXP p){
    covafill<double>* cfp = new covafill<double>(asMatrix(coord),
						       asVector(obs),
						       asVector(h),
						       asInteger(p));

    if(cfp == NULL){
      return R_NilValue;
    }
    SEXP val = R_MakeExternalPtr(cfp, install("covafillPointer"), R_NilValue);
    return val;
  }

  SEXP getFillDim(SEXP sp){
    if(R_ExternalPtrTag(sp) != install("covafillPointer"))
      Rf_error("The pointer must be to a covafill object");   
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);
    return asSEXP(ptr->getDim());
  }

  SEXP getFillDegree(SEXP sp){
    if(R_ExternalPtrTag(sp) != install("covafillPointer"))
      Rf_error("The pointer must be to a covafill object");   
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);
    return asSEXP(ptr->p);
  }

  SEXP getFillBandwith(SEXP sp){
    if(R_ExternalPtrTag(sp) != install("covafillPointer"))
      Rf_error("The pointer must be to a covafill object");   
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);
    return asSEXP(ptr->h);
  }

  SEXP setFillBandwith(SEXP sp, SEXP h){
    
    if(R_ExternalPtrTag(sp) != install("covafillPointer"))
      Rf_error("The pointer must be to a covafill object");   
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);
    if(LENGTH(h) == 1){
      ptr->setH(asDouble(h));
    }else{
      ptr->setH(asVector(h));
    }
    int res = 1;
    return asSEXP(res);
  }

  SEXP predictFill(SEXP sp, SEXP x){
    if(R_ExternalPtrTag(sp) != install("covafillPointer"))
      Rf_error("The pointer must be to a covafill object");   
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);

    if(isMatrix(x)){
      int lsdim = 1 + ptr->getDim();
      if(ptr->p >= 2)
	lsdim += 0.5 * ptr->getDim() * (ptr->getDim() + 1);
      if(ptr->p >= 3)
	lsdim += (ptr->p - 2) * ptr->getDim();
      MatrixXd res(nrows(x),lsdim);
      MatrixXd x0 = asMatrix(x);
      for(int i = 0; i < nrows(x); ++i)
	res.row(i) = ptr->operator()((vector)x0.row(i), true);
      return asSEXP(res);
    }else if(isNumeric(x)){
      return asSEXP(ptr->operator()(asVector(x), true));
    }else{
      error("Element must be a matrix or numeric vector");
    }
    return R_NilValue;
  }


  SEXP lnoResiduals(SEXP sp, SEXP excludeRadius){
    if(!(isNumeric(excludeRadius) && LENGTH(excludeRadius) == 1))
      Rf_error("Exclude radius must be a scalar");
    if(R_ExternalPtrTag(sp) != install("covafillPointer"))
      Rf_error("The pointer must be to a covafill object");
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);

    double er = asDouble(excludeRadius);
    MatrixXd x0 = ptr->coordinates;
    vector y0 = ptr->observations;
    vector res(y0.size());

    for(int i = 0; i < x0.rows(); ++i)
      res.row(i) = ptr->operator()((vector)x0.row(i),er) - y0(i);

    return asSEXP(res);
  }
}
