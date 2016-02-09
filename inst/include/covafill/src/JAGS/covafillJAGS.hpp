

namespace jags {
  namespace covafillJAGS {

    class covafillJAGS : public ArrayFunction
    {
    public:
	covafillJAGS();
	void evaluate(double *value, std::vector<double const *> const &args,
		      std::vector<std::vector<unsigned int> > const &dims) 
	    const;
	std::vector<unsigned int> 
	    dim(std::vector<std::vector<unsigned int> > const &dims,
		std::vector<double const *> const &values) const;
	bool checkParameterDim(std::vector <std::vector<unsigned int> > const &dims) const;

    };

  }}
