#ifndef MULTIPOISSONGAMMAMODEL_H
#define MULTIPOISSONGAMMAMODEL_H
//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// See: http://arxiv.org/abs/arXiv:1002.1111
// 
// File: MultiPoissonGammaModel.h
// Description: Implement the Poisson-Gamma model 
//              marginalized over the nuisance parameters).
// 
// Created: 11-Jun-2010
// Updated: 04-Jul-2015 HBP Renamed MultiPoissonGammaModel.cc 
//
//--------------------------------------------------------------
#include <vector>
#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"
#include "PDFunction.h"

/** Implement the Poisson-Gamma model marginalized 
    over the nuisance parameters.
     <p>
        The Poisson mean (per bin) is defined by,   
        \f[
	\sigma \, \epsilon + \mu,
	\f]
	where \f$\sigma\f$ is the signal cross section, \f$\epsilon\f$
	the effective luminosity (that is, the signal efficiency times 
	acceptance times the integrated luminosity), and \f$\mu\f$ 
	the background.  
        <p>
	The priors for \f$\epsilon\f$ and \f$\mu\f$ are modeled
	as gamma densities,
        \f[
	f(z, \gamma, \beta) = \frac{(z /\beta)^{\gamma - 1} e^{-z/\beta}}
	{\Gamma(\gamma) \beta},
	\f]
	where, for the background, 
	\f{eqnarray*}{
	\gamma & = & y + 1/2, \\
	\beta  & = & 1 / b,
	\f}
	and for the signal, 
	\f{eqnarray*}{
	\gamma & = & x + 1/2, \\
	\beta  & = & 1 / a.
	\f}
	The counts \f$y\f$ and \f$x\f$ are defined so that 
      \f{eqnarray*}{
      E[y] & = & b \, \mu, \\ 
      E[x] & = & a \, \epsilon,
      \f}
      where \f$a\f$ and \f$b\f$ are scale factors.
      This prior can be motivated as follows. Given a count \f$y\f$
      and scale factor \f$b\f$, the likelihood and associated reference prior
      for \f$\mu\f$ are
      \f{eqnarray*}{
      p(y | b \, \mu) & = & \textrm{Poisson}(y, b \, \mu), \\
                   & = & \frac{(b \, \mu)^y e^{-b \, \mu}}{\Gamma(y + 1)},
      \quad\textrm{and}\\			      \
      \pi(\mu) & = & 1 / \sqrt{\mu},
      \f}
      which, via Bayes theorem, yields a gamma posterior density for \f$\mu\f$.
      The latter
      serves as the prior for
      \f$\mu\f$ in the subsequent inference. An analogous consideration 
      applies to the effective luminosity prior.
      <p>
      Often the counts \f$y\f$ and \f$x\f$ and associated
      scale factors are not given. For example, for the background
      one may have instead 
      estimates in the form \f$c \pm \delta c\f$ for the background
      \f$\mu\f$. One can 
      infer an effective count \f$y\f$ and a scale factor \f$b\f$ by
      identifying 
      \f{eqnarray*}{
      c & = & \beta \, (\gamma - 1) \quad\textrm{as the mode, and}\\
      \delta c^2 & = & \beta^2 \, \gamma \quad\textrm{as the variance}
      \f}
      of a gamma density
	where the count \f$y\f$ and scale factor \f$b\f$ are given by 
	\f{eqnarray*}{
	y + 1 & \equiv & \gamma, \quad\textrm{and} \\
	b & \equiv & 1 / \beta,
	\f}
	with 
      \f{eqnarray*}{
      \gamma & = & [k+2 + \sqrt{(k+2)^2 - 4}]/2, \\
      \beta  & = & [\sqrt{c^2 + 4 \, \delta c^2} - c] / 2, 
      \f}
      where \f$k = (c / \delta c)^2\f$. And similarly for the effective
      luminosities.
*/
class MultiPoissonGammaModel : public PDFunction
{
 public:
  
  ///
  MultiPoissonGammaModel();
  
  /** Implement marginalized Poisson-Gamma model.
      We use the notation of the paper http://arxiv.org/abs/arXiv:1002.1111
      @param data - observed counts
      @param x - counts for effective luminosity prior
      @param a - scale factors for effective luminosity prior
      @param y - counts for background prior
      @param b - scale factors for background prior
      @maxcount - maximum allowed value of observed count 
      (this sets the sizes of some internal arrays).
  */
  MultiPoissonGammaModel(std::vector<double>& data,
			 std::vector<double>& x, 
			 std::vector<double>& a,
			 std::vector<double>& y, 
			 std::vector<double>& b,
			 int maxcount=100000);
  
  /** Implement marginalized Poisson-Gamma model.
      @param data - observed counts
      @param x - counts for effective luminosity prior
      @param a - scale factor for effective luminosity prior
      @param y - counts for background prior
      @param b - scale factor for background prior
      @maxcount - maximum allowed value of observed count 
      (this sets the sizes of some internal arrays).
  */
  MultiPoissonGammaModel(std::vector<double>& data,
			 std::vector<double>& x,
			 double a,
			 std::vector<double>& y, 
			 double b,  
			 int maxcount=100000);

  /** Implement marginalized Poisson-Gamma model.
      @param data - observed count
      @param x - count for effective luminosity prior
      @param a - scale factor for effective luminosity prior
      @param y - count for background prior
      @param b - scale factor for background prior
      @maxcount - maximum allowed value of observed count 
      (this sets the sizes of some internal arrays).
  */
  MultiPoissonGammaModel(double data,
			 double x,
			 double a,
			 double y, 
			 double b,  
			 int maxcount=100000); 
 
  ///
  ~MultiPoissonGammaModel();
  
  ///
  void setX(std::vector<double>& x, std::vector<double>& a)
  {
    _x = x;
    _a = a;
  }
  
  ///
  void setY(std::vector<double>& y, std::vector<double>& b)
  {
    _y = y;
    _b = b;
  }    
  
  /** Generate data for one experiment.
        @param sigma - value of parameter of interest
  */
  std::vector<double>& generate(double sigma);
  
  /** Compute likelihood.
      @param data  - observed counts
      @param sigma - parameter of interest 
  */
  double operator() (std::vector<double>& data, double sigma);
  
  /** Compute likelihood using cached data.
      @param sigma - parameter of interest 
  */
  double operator() (double sigma);
    
  std::vector<double>& counts() { return _data; }
  
 private:
  std::vector<double> _data;
  std::vector<double> _x;
  std::vector<double> _a;
  std::vector<double> _y;
  std::vector<double> _b;
  int _maxcount;
  ROOT::Math::Random<ROOT::Math::GSLRngMT> _gslRan;
};

#endif

