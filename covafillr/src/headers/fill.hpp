
#define CHECK_FILL_POINTER(name)		\
  if(R_ExternalPtrTag(name) != Rf_install("covafillPointer")) \
    Rf_error("The pointer must be to a covafill object"); \
  if(!R_ExternalPtrAddr(name)) \
    Rf_error("The pointer address is not valid");

extern "C" {

  static void finalizeFill(SEXP ptr){
    if(!R_ExternalPtrAddr(ptr))
      return;
    delete (covafill<double>*)R_ExternalPtrAddr(ptr);
    R_ClearExternalPtr(ptr);    
  }
  
  SEXP MakeFill(SEXP coord,SEXP obs,SEXP h,SEXP p){
    covafill<double>* cfp = new covafill<double>(asMatrix(coord),
						       asVector(obs),
						       asVector(h),
						       asInteger(p));

    if(cfp == NULL){
      return R_NilValue;
    }
    SEXP val = R_MakeExternalPtr(cfp, Rf_install("covafillPointer"), R_NilValue);
    PROTECT(val);
    R_RegisterCFinalizerEx(val, finalizeFill, TRUE);
    UNPROTECT(1);
    return val;
  }

  SEXP getFillDim(SEXP sp){
    CHECK_FILL_POINTER(sp);
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);
    return asSEXP(ptr->getDim());
  }

  SEXP getFillDegree(SEXP sp){
    CHECK_FILL_POINTER(sp);
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);
    return asSEXP(ptr->p);
  }

  SEXP getFillBandwith(SEXP sp){
    CHECK_FILL_POINTER(sp);
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);
    return asSEXP(ptr->h);
  }

  SEXP getFillObservations(SEXP sp){
    CHECK_FILL_POINTER(sp);
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);
    return asSEXP(ptr->observations);
  }

  SEXP getFillCoordinates(SEXP sp){
    CHECK_FILL_POINTER(sp);
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);
    return asSEXP(ptr->coordinates);
  }

  SEXP setFillBandwith(SEXP sp, SEXP h){
    CHECK_FILL_POINTER(sp);
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
    CHECK_FILL_POINTER(sp);
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);

    if(Rf_isMatrix(x)){
      int lsdim = 1;
      if(ptr->p >= 1)
	lsdim += ptr->getDim();
      if(ptr->p >= 2)
	lsdim += 0.5 * ptr->getDim() * (ptr->getDim() + 1);
      if(ptr->p >= 3)
	lsdim += (ptr->p - 2) * ptr->getDim();
      MatrixXd res(Rf_nrows(x),lsdim);
      MatrixXd x0 = asMatrix(x);
      for(int i = 0; i < Rf_nrows(x); ++i)
	res.row(i) = ptr->operator()((vector)x0.row(i), true);
      return asSEXP(res);
    }else if(Rf_isNumeric(x)){
      return asSEXP(ptr->operator()(asVector(x), true));
    }else{
      Rf_error("Element must be a matrix or numeric vector");
    }
    return R_NilValue;
  }


  SEXP predictFillSE(SEXP sp, SEXP x){
    CHECK_FILL_POINTER(sp);   
    covafill<double>* ptr=(covafill<double>*)R_ExternalPtrAddr(sp);
    
    if(Rf_isMatrix(x)){
      MatrixXd x0 = asMatrix(x);
      
      int lsdim = 1;
      if(ptr->p >= 1)
	lsdim += ptr->getDim();
      if(ptr->p >= 2)
	lsdim += 0.5 * ptr->getDim() * (ptr->getDim() + 1);
      if(ptr->p >= 3)
	lsdim += (ptr->p - 2) * ptr->getDim();

      MatrixXd res(Rf_nrows(x),lsdim);
      MatrixXd resSE(Rf_nrows(x),lsdim);
      
      Array<Array<double,Dynamic,1>, Dynamic,1> tmp(2);
      for(int i = 0; i < Rf_nrows(x); ++i){
	tmp = ptr->operator()((vector)x0.row(i),int(0), true);
	res.row(i) = tmp(0);
	resSE.row(i) = tmp(1);
      }

      SEXP vecOut = PROTECT(Rf_allocVector(VECSXP, 2));
      SEXP sr1 = PROTECT(asSEXP(res));
      SEXP sr2 = PROTECT(asSEXP(resSE));
      SET_VECTOR_ELT(vecOut,0,sr1);
      SET_VECTOR_ELT(vecOut,1,sr2);

      UNPROTECT(3);
      return vecOut;
      
    }else{
      Rf_error("Element must be a matrix or numeric vector");
    }
    return R_NilValue;
  }



  SEXP lnoResiduals(SEXP sp, SEXP excludeRadius){
    CHECK_FILL_POINTER(sp);
    if(!(Rf_isNumeric(excludeRadius) && Rf_length(excludeRadius) == 1))
      Rf_error("Exclude radius must be a scalar");
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
