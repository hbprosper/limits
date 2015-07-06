#ifndef MULTIPOISSONGAMMA_H
#define MULTIPOISSONGAMMA_H
//--------------------------------------------------------------
// File: MultiPoissonGamma.h
// Description: Implement the multi-Poisson-gamma model averaged
//              with respect to a prior specified as a swarm
//              of points. Model assumes mean = eff*L + bkg for
//              each bin.
// 
// Created: 11-Jun-2010 Harrison B. Prosper & Supriya Jain
// Updated: 11-Aug-2014 HBP - add option to profile model
//                      model (not yet implemented!).
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
   p(\sigma | D) & = & \int \textrm{PoissonGamma}(D | \sigma, y, x, b, a) 
   \, \pi(y, x, b, a) \, dy \, dx \, db \, da, \\
   & \approx & \frac{1}{K}
   \sum_{i=1}^K \textrm{PoissonGamma}(D | \sigma, y_i, x_i, b_i, a_i), 
   \f}
   where \f$\pi(y, x, b, a)\f$ is an evidenced-based prior specified
   as a swarm of \f$K\f$ sampled points \f$\{(y_i, x_i, b_i, a_i)\}\f$. The
   parameters \f$y\f$ and \f$x\f$ are backgrounds and effective luminosities
   (signal efficiency times acceptance times luminosities), respectively,
   with associated scale factors  \f$b\f$, and \f$a\f$. 
*/
class MultiPoissonGamma : public PDFunction
{
 public:
  ///
  MultiPoissonGamma();

  /** Default constructor. 
      @param filename - name of text file containing counts, etc.
  */
  MultiPoissonGamma(std::string filename);
  
  /** Main constructor.  
      @param N - observed counts
  */
  MultiPoissonGamma(std::vector<double>& N);
 
  ///
  ~MultiPoissonGamma();
	     
  /** Generate data for one experiment.
      @param sigma - value of parameter of interest
  */
  std::vector<double>& generate(double sigma);
  
  /** Compute likelihood.
      @param N - observed data
      @param sigma - value of parameter of interest 
  */
  double operator() (std::vector<double>& N, double sigma);

  /** Compute likelihood using internally cached data.
      @param sigma - value of parameter of interest 
  */
  double operator() (double sigma);

  /** If true, profile rather than average.
      Not yet implemented.
   */
  void profile(bool yes=true) {_profile=yes;}

  /** Add one set of backgrounds and effective luminosities.
      @param bkg  - backgrounds
      @param dbkg - background uncertainties
      @param efl  - effective luminosities
      @param defl - effective luminosity uncertainties

     <p>
     This class uses MultiPoissonGammaModel, which requires
     counts \f$y\f$ and \f$x\f$ and their associated
      scale factors \f$b\f$ and \f$a\f$, respectively, that define
      the gamma priors for the backgrounds and effective luminosities.
      However, often the evidence-based prior is specified as
      estimates in the form \f$c \pm \delta c\f$ for the effective
      luminosities and backgrounds. The add method converts these 
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
      applies to the effective luminosities.
   */
  void add(std::vector<double>& bkg, std::vector<double>& dbkg,
	   std::vector<double>& efl, std::vector<double>& defl);
	   
  ///
  void update(int ii, std::vector<double>& eff, std::vector<double>& deff);

  ///
  void set(int ii);
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
    int _nbins;
    int _index;
    bool _profile;

    void _convert(std::vector<double>& efl, std::vector<double>& defl,
		  std::vector<double>& x,   std::vector<double>& a);    
};

#endif

