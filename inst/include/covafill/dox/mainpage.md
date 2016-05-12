# covafill: Local Polynomial Regression of State Dependent Covariates in State-Space Models

covafill is a C++ template library for local polynomial regression of covariates in state-space models. The covafill library is based on the [Eigen](http://http://eigen.tuxfamily.org) library for linear algebra, and includes several modules:

- The [Core module](@ref core) which provides the base functionality for local polynomial regression
- The [Tree module](@ref tree) which provides a search tree approximation to local polynomial regression
- The [Interpolate module](@ref interpolate) which provides classes for cubic interpolation in 1-3 dimensions. This module is only inteded for internal use.
- The [JAGS module](@ref jags), which provides a module for using covafill with [JAGS](http://mcmc-jags.sourceforge.net/)
- The [TMB module](@ref tmb) which provides functionality to use covafill with [TMB](http://tmb-project.org).

## The Core module

The Core module provides the class [covafill](@ref covafill) for local polynomial regression.

### Local polynomial regression

For simplicity, consider the univariate model

\f[
y_i = g(x_i) + \epsilon_i
\f]

where \f$g:\mathbb{R}\mapsto\mathbb{R}\f$ is a smooth function and \f$ \epsilon_i\sim N(0,\sigma^2)\f$.
To do local polynomial regression of \f$g\f$ at \f$x_0\f$, we do a Taylor expansion of order \f$p\f$,

\f[
g(x) \approx g(x_0) + g^{(1)}(x_0)(x-x_0)  + \frac{1}{2!} g^{(2)}(x_0)(x-x_0)^2 + \cdots + \frac{1}{p!} g^{(p)}(x_0)(x-x_0)^p
\f]
Substituting into the original model,
\f[
y_i = g(x_0) + g^{(1)}(x_0)(x-x_0)  + \frac{1}{2!} g^{(2)}(x_0)(x-x_0)^2 + \cdots + \frac{1}{p!} g^{(p)}(x_0)(x-x_0)^p + \epsilon_i
\f]
we obtain a linear model with coefficients \f$ \mathbf{\theta} = (g(x_0), g^{(1)}(x_0), g^{(2)}(x_0), \ldots, g^{(p)}(x_0))^T \f$, obervations \f$ \mathbf{Y} = (y_1, y_2, \ldots, y_n)^T \f$, and the design matrix
\f[
\mathbf{X} = \left(\begin{array}{ccccc}
1 & (x_1-x_0) & \frac{1}{2!} (x_1-x_0)^2 & \cdots & \frac{1}{p!} (x_1-x_0)^p \\
1 & (x_2-x_0) & \frac{1}{2!} (x_2-x_0)^2 & \cdots & \frac{1}{p!} (x_2-x_0)^p \\
\vdots & \vdots & \vdots &   & \vdots \\
1 & (x_n-x_0) & \frac{1}{2!} (x_n-x_0)^2 & \cdots & \frac{1}{p!} (x_n-x_0)^p
\end{array}\right)
\f]
As we are interested in a local estimate, observations are weighed by their distance to \f$x_0\f$. The weights form the diagonal matrix \f$\mathbf{W}\f$ with
\f[
w_{ii} = \det(H^{-1}) \left(1 - \| H^{-1} \cdot ( x_i - x_0) \| ^ 2 \right) \vee 0
\f]

Now the estimates are obtained by
\f[
\hat{\mathbf{\theta}} = (\mathbf{X}^T\mathbf{W}\mathbf{X})^{-1} \mathbf{X}^T\mathbf{W}\mathbf{Y}
\f]
giving both the estimated function value at \f$ x_0 \f$ and estimates of the first \f$ p \f$ derivatives.

A covariance matrix for the estimates can be calculated by

\f[
V(\hat{\mathbf{\theta}}) = (\mathbf{X}^T\mathbf{W}\mathbf{X})^{-1} (\mathbf{R}^T \mathbf{W} \mathbf{R}) (N-q)^{-1}
\f]

where \f$ N \f$ is the number of observations with non-negative weights, \f$ q \f$ is the number of regressors, and

\f[
\mathbf{R} = \mathbf{Y}-\mathbf{X}\hat{\mathbf{\theta}}
\f] 


Note that it is not necessary to have a properly normalized kernel function as the normalizing constant vanishes in the calculation of both the estimates and covariance matrix.


## The Tree module

The Tree module contains the covatree class for a search tree approximation to local polynomial regression.
The covatree class builds a simple search tree by splitting the data in the coordinate with the highest variance.
The split is performed at the mean of the coordinates. 
At terminal notes, the class does local polynomial regression at the corners of the bounding box of the coordinates related to the note and calculates the necessary coefficients to do cubic interpolation between the corners.


## The JAGS and TMB modules

The JAGS and TMB modules provide functionality to use covafill within these tools.
The JAGS module provides a JAGS Module including a function, covafill, to call in a JAGS model (See example below).
The TMB module include functions to evaluate a covafill or covatree object from a TMB model such that the estimated derivatives are used in the automatic differentiation (See example below).


### JAGS example

~~~~~~~~~~~~~{.cpp}
model {
      cf <- covafill(x,obsC,obs,h,2.0)
      sigma ~ dunif(0,100)
      tau <- pow(sigma, -2)
      for(i in 1:N) {
      	    y[i] ~ dnorm(cf[i],tau)
	    }
}
~~~~~~~~~~~~~

### TMB example

~~~~~~~~~~~~~{.cpp}
#include <TMB.hpp>
#include <covafill/TMB>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(obs);
  DATA_MATRIX(coord);
  DATA_VECTOR(covObs);
  DATA_INTEGER(p);
  DATA_VECTOR(h);

  PARAMETER(logObsSd);
  PARAMETER(logObsTSd);
  PARAMETER(logStatSd);
  PARAMETER_MATRIX(x);

  Type nll = 0.0;
  covafill<Type> cf(coord,covObs,h,p);

  // Contribution from states
  for(int i = 1; i < x.cols(); ++i){
    nll -= dnorm(x(0,i), x(0,i-1), exp(logStatSd),true);
    nll -= dnorm(x(1,i), x(1,i-1), exp(logStatSd),true);
  }

  // contribution from observations
  for(int i = 0; i < obs.cols(); ++i){
    nll -= dnorm(obs(0,i), x(0,i), exp(logObsSd),true);
    nll -= dnorm(obs(1,i), x(1,i), exp(logObsSd),true);
    vector<Type> tmp = x.col(i);
    Type val = evalFill((CppAD::vector<Type>)tmp, cf)[0];
    nll -= dnorm(obs(2,i), val, exp(logObsTSd),true);
  }
 
    
  return nll;
}
~~~~~~~~~~~~~
