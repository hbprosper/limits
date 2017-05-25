#ifndef Bayes_H
#define Bayes_H
//--------------------------------------------------------------
//
// File: Bayes.cc
// Description: Standalone Bayes limit calculator.
// 
// Created: 11 Jan 2011 Harrison B. Prosper
// Updated: 06 Mar 2014 HBP clean up
//          30 May 2015 HBP add a constructor to take RooFit
//                      models directly rather than through
//                      PDFWrapper
//
//--------------------------------------------------------------
#include <vector>
#include <string>
#include "PDFunction.h"
#include "PriorFunction.h"
#include "Math/Interpolator.h"
#ifdef __WITH_ROOFIT__
#include "PDFWrapper.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#endif
//--------------------------------------------------------------
/** Compute Bayesian limits.
 */
class Bayes
{
public:
  Bayes () {}

  /** Compute Bayesian limits.
      @param model  - probability density function (pdf)
      @param data   - observed data
      @param poimin - minimum of parameter of interest
      @param poimax - maximum of parameter of interest
      @param cl     - confidence level
      @param prior  - prior funtion (default = flat)
  */
  Bayes(PDFunction& model,
	std::vector<double>& data,
	double poimin,
	double poimax,
	double cl=0.95,
	PriorFunction* prior_=0);
  
#ifdef __WITH_ROOFIT__
    /** RooFit constructor.  
      @param pdf    - pdf object of type RooAbsPdf
      @param obs    - set of data objects
      @param poi    - parameter of interest
      @param cl     - confidence level
      @param prior  - prior function of type RooAbsPdf (default = flat)
   */
  Bayes(RooAbsPdf& pdf, RooArgSet& obs, RooRealVar& poi,
	double cl=0.95,
	RooAbsPdf* prior_=0);
#endif
  
  virtual ~Bayes();

  /** Compute prior.
   */  
  double prior(double poi);
  
  /** Compute likelihood.
   */  
  double likelihood(double poi);

  /** Compute normalization of posterior density.
   */  
  double normalize();

  /** Compute posterior density at specified value of parameter of interest.
   */  
  double posterior(double poi);

  /// Compute cdf of posterior density.
  double cdf(double poi);

  /// Compute p-quantile of posterior density.
  double quantile(double p=-1);

  /// Compute MAximum Posterior estimate (mode of posterior density).
  void MAP(double CL, std::vector<double>& results);
  
  void setData(std::vector<double>& d) { _data = d; normalize(); }
  
  std::vector<double>& data() {return _data;}

  double CL() { return _cl; }
  
  /// Compute Bayesian Z=sign(B10)*sqrt(2*|B10|), where B10 is the Bayes factor.
  double zvalue(double mu=1);

private:
  PDFunction*    _pdf;
  std::vector<double> _data;
  double _poimin;
  double _poimax;
  double _cl;
  PriorFunction* _prior;
#ifdef __WITH_ROOFIT__
  RooAbsPdf*     _rfprior;
  RooRealVar*    _rfpoi;
#endif
  bool   _normalize;
  
  ROOT::Math::Interpolator* _interp;  
  int _nsteps;
  std::vector<double> _x;
  std::vector<double> _y;
  
  double _normalization;
  double _likeprior(double poi);
  double _q(double prob);
  double _f(double prob);
  double _nsig;
  double _postmax;
};

#endif


