/** \defgroup jags JAGS module
*
* The JAGS module of covafill provides the implementation of a JAGS user module with the function covafill for local polynomial regression in a JAGS model.
* \verbatim
#include <covafill/JAGS>
\endverbatim
*/

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
