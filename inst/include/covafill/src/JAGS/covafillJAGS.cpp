
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
