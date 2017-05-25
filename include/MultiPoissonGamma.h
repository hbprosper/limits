#ifndef MULTIPOISSONGAMMA_H
#define MULTIPOISSONGAMMA_H
//--------------------------------------------------------------
// File: MultiPoissonGamma.h
// Description: Implement the multi-Poisson-gamma model averaged
//              with respect to a prior specified as a swarm
//              of points. Model assumes mean = mu*s + b for
//              each bin, where mu is the signal strength.
// 
// Created: 11-Jun-2010 Harrison B. Prosper & Supriya Jain
// Updated: 11-Aug-2014 HBP - add option to profile model
//                      model (not yet implemented!).
//          21-Jun-2016 HBP - add histogram reading option
//          25-May-2017 HBP - use S and B instead of efl and bkg!
//--------------------------------------------------------------
#include <vector>
#include <algorithm>
#include "TRandom3.h"
#include "PDFunction.h"
#include "MultiPoissonGammaModel.h"

/** Implement the MultiPoissonGammaModel averaged over an evidence-based prior.
   <p>
   The class computes
   \f{eqnarray*}{
   p(\mu | D) & = & \int \textrm{PoissonGamma}(D | \mu, x, y, a, b) 
   \, \pi(x, y, a, b) \, dx \, dy \, da \, db, \\
   & \approx & \frac{1}{K}
   \sum_{i=1}^K \textrm{PoissonGamma}(D | \mu, x_i, y_i, a_i, b_i), 
   \f}
   where \f$\pi(x, y, a, b)\f$ is an evidenced-based prior specified
   as a swarm of \f$K\f$ sampled points \f$\{(x_i, y_i, a_i, b_i)\}\f$. The
   parameters \f$x\f$ and \f$y\f$ are the signals and backgrounds, respectively,
   with associated scale factors  \f$a\f$, and \f$b\f$. 
*/
class MultiPoissonGamma : public PDFunction
{
 public:
  ///
  MultiPoissonGamma();

  /** 
      @param filename - name of text file containing counts, or histogram names.
  */
  MultiPoissonGamma(std::string filename);
  
  /** 
      @Param N - observed counts
  */
  MultiPoissonGamma(std::vector<double>& N);
  
  ///
  virtual ~MultiPoissonGamma();
	     
  /** Generate data for one experiment.
      @param mu - signal strength (parameter of interest)
  */
  std::vector<double>& generate(double mu);
  
  /** Compute likelihood.
      @param N - observed data
      @param mu - signal strength (parameter of interest)
  */
  double operator() (std::vector<double>& N, double mu);

  /** Compute likelihood using internally cached data.
      @param mu - signal strength (parameter of interest) 
  */
  double operator() (double mu);

  /** If true, profile rather than average.
      Not yet implemented.
   */
  void profile(bool yes=true) {_profile=yes;}

  /** Add one set of signals and backgrounds and associated uncertainties.
      @param sig  - signals
      @param dsig - signal uncertainties
      @param bkg  - backgrounds
      @param dbkg - background uncertainties

     <p>
     This class uses MultiPoissonGammaModel, which requires
     counts \f$x\f$ and \f$y\f$ and their associated
      scale factors \f$a\f$ and \f$b\f$, respectively, that define
      the gamma priors for the backgrounds and signals.
      However, often the evidence-based prior is specified as
      estimates in the form \f$c \pm \delta c\f$ for the background
      and signals. The add method converts these 
      estimates into effective counts and scale factors and then
      instantiates a MultiPoissonGammaModel.
      <p>
      Consider the backgrounds.
      An effective count \f$y\f$ and scale factor \f$b\f$ is found by
      identifying 
      \f{eqnarray*}{
      c & = & \beta \, (\gamma - 1) 
      \quad\textrm{as the mode, and}\\
      \delta c^2 & = & \beta^2 \, \gamma \quad\textrm{as the variance}
      \f}
      of the gamma density
       \f[
	f(\mu, \gamma, \beta) = \frac{(\mu / \beta)^{\gamma - 1} e^{-\mu/\beta}}
	{\Gamma(\gamma) \beta },
	\f]
	where 
	\f{eqnarray*}{
	y + 1 & \equiv & \gamma, \quad\textrm{and} \\
	b & \equiv & 1 / \beta,
	\f}
       with,
      \f{eqnarray*}{
      \gamma & = & [k+2 + \sqrt{(k+2)^2 - 4}]/2, \\
      \beta  & = & [\sqrt{c^2 + 4 \, \delta c^2} - c] / 2, 
      \f}
      where \f$k = (c / \delta c)^2\f$. An analogous consideration
      applies to the signals.
   */
  void add(std::vector<double>& sig, std::vector<double>& dsig,
	   std::vector<double>& bkg, std::vector<double>& dbkg);
	   
  ///
  void update(int ii, std::vector<double>& sig, std::vector<double>& dsig);

  ///
  void set(int ii);

  ///
  void reset();

  ///
  void setSeed(int seed);

  ///
  std::vector<double> counts() { return _N; }

  /// Sample size.
  int size() { return _model.size(); }
  
 private:
    std::vector<double> _N;
    std::vector<double> _Ngen;
    std::vector<MultiPoissonGammaModel> _model;
    
    TRandom3 _random;
    int  _nbins;
    int  _index;
    bool _profile;

    void _convert(std::vector<double>& sig, std::vector<double>& dsig,
		  std::vector<double>& x,   std::vector<double>& a);
    void _readTextFile(std::vector<std::string>& records);
    void _readRootFile(std::vector<std::string>& records);
    void _getcounts(std::string rfilename, std::string histname,
		    std::vector<double>& c, std::vector<double>& dc);    
};

#endif

