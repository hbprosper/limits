#ifndef MULTIPOISSONGAMMA_H
#define MULTIPOISSONGAMMA_H
//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: MultiPoissonGamma.h
// Description: Implement the Poisson-Gamma model 
//              marginalized over the nuisance parameters).
// 
// Created: 11-Jun-2010
// Updated: 04-Jul-2015 HBP Renamed MultiPoissonGamma.cc 
//
//--------------------------------------------------------------
#include <vector>
#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"
#include "PDFunction.h"

/** Implement the Poisson-Gamma model marginalized 
    over the nuisance parameters.
*/
class MultiPoissonGamma : public PDFunction
{
 public:
  
  ///
  MultiPoissonGamma();
  
  /** Constructor:  
      @param yy - counts for background prior
      @param xx - counts for signal prior
      @param b  - scale for background prior
      @param a  - scale for signal prior
      @maxcount - maximum allowed value of observed count 
      (this sets the sizes of some internal arrays).
      <p>
        The Poisson mean is defined by:   
        \f$$\sigma*\epsilon + b.\f$$
        <p>
	The priors for \f$\epsilon\f$ and \f$\mu\f$ are modeled
	as gamma densities,
        \f$$f(x, \gamma, \beta) = \frac{(x /\beta)^{\gamma - 1} e^{-x/\beta}}
	{\Gamma(\gamma) \beta},\f$$
	where, for the background, \f$\gamma = yy + 1/2\f$ and \f$\beta = 1/b\f$,
	and for the signal, \f$\gamma = xx + 1/2\f$ and \f$\beta = 1/a\f$.    
    */
    MultiPoissonGamma(std::vector<double>& yy, 
		      std::vector<double>& xx,
		      std::vector<double>& b, 
		      std::vector<double>& a,
		      int maxcount=100000);
    /** Constructor:  
        @param yy - counts for background prior
        @param xx - counts for signal prior
        @param b  - scale for background prior
        @param a  - scale for signal prior
        @maxcount - maximum allowed value of observed count 
    */
    MultiPoissonGamma(std::vector<double>& yy,  
		      std::vector<double>& xx, 
		      double b,  
		      double a,
		      int maxcount=100000); 
  
    ///
    ~MultiPoissonGamma();
    
    /** Generate data for one experiment.
        @param sigma - value of parameter of interest
    */
    std::vector<double>& generate(double sigma);
	     
    /** Compute PDF.
        @param xdata - observed counts
        @param sigma - parameter of interest 
    */
    double operator() (std::vector<double>& xdata, double sigma);
    
 private:	     
    std::vector<double> _yy;
    std::vector<double> _xx;

    int _nbins;
    
    std::vector<double> _b;
    std::vector<double> _a;
    int _maxcount;
    std::vector<double> _xdata;
    ROOT::Math::Random<ROOT::Math::GSLRngMT> _gslRan;
};

#endif

