#ifndef Bayes_H
#define Bayes_H
//--------------------------------------------------------------
//
// File: Bayes.cc
// Description: Standalone Bayes limit calculator.
// 
// Created: 11 Jan 2011 Harrison B. Prosper
// Updated: 06 Mar 2014 HBP clean up
//          11 Feb 2016 HBP rename shadowed variable
//--------------------------------------------------------------
#include <vector>
#include <string>
#include "PDFunction.h"
#include "PriorFunction.h"
#include "Math/Interpolator.h"
//--------------------------------------------------------------
class Bayes
{
public:
  Bayes () {}

  /** 
      @param model  - pdf
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
	double cl_=0.95,
	PriorFunction* prior_=0);

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

private:
  PDFunction*    _pdf;
  std::vector<double> _data;
  double _poimin;
  double _poimax;
  double _cl;
  PriorFunction* _prior;
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


