//--------------------------------------------------------------
// File: Bayes.cc
// Description: Implement standalone Bayes calculator
// 
// Created: 11 Jan 2011 HBP
// Updated: 13 Mar 2011 HBP - fix normalization of post.
//          06 Mar 2014 HBP
//          30 May 2015 HBP - implement direct RooFit interface.
//          03 May 2018 HBP - add estimate and uncertainty methods
//--------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <stdlib.h>

#include "Bayes.h"
#include "TMinuit.h"
#include "TMath.h"
#include "Math/WrappedFunction.h"
#include "Math/Integrator.h"
#include "Math/RootFinder.h"

using namespace std;
// ---------------------------------------------------------------------------

// function to be minimized
namespace {
  const int MAXITER=10000;
  const double TOLERANCE=1.e-5;
  Bayes* OBJ=0;
  void nlpFunc(int&    /*npar*/, 
	       double* /*grad*/, 
	       double& fval,
	       double* xval,
	       int     /*iflag*/)
  {
    double poi = xval[0];
    fval = -log(OBJ->posterior(poi));
  }
};

Bayes::Bayes(PDFunction& model,
	     std::vector<double>& d,
	     double poimin,
	     double poimax,
	     double cl_,
	     PriorFunction* prior_)
  : _pdf(&model),
    _data(d),
    _poimin(poimin),
    _poimax(poimax),
    _cl(cl_),
    _prior(prior_),
#ifdef __WITH_ROOFIT__
    _rfprior(0),
    _rfpoi(0),
#endif
    _normalize(true),
    _interp(0),
    _nsteps(50),
    _x(vector<double>()),
    _y(vector<double>()),
    _MAPdone(false),
    _result(std::pair<double, double>(0, 0)),
    _verbosity(-1)
{
  if ( getenv("limits_verbosity") != (char*)0 )
    _verbosity = atoi(getenv("limits_verbosity"));

  assert( _poimax > _poimin ); 
  OBJ = this;
  normalize();
}

#ifdef __WITH_ROOFIT__
Bayes::Bayes(RooAbsPdf& pdf, RooArgSet& obs, RooRealVar& poi,
	     double cl_,
	     RooAbsPdf* prior_)
  : _pdf(new PDFWrapper(pdf, obs, poi)),
    _data(vector<double>(obs.getSize())),
    _poimin(poi.getMin()),
    _poimax(poi.getMax()),
    _cl(cl_),
    _prior(0),
    _rfprior(prior_),
    _rfpoi(&poi),
    _normalize(true),
    _interp(0),
    _nsteps(50),
    _x(vector<double>()),
    _y(vector<double>()),
    _MAPdone(false),
    _result(std::pair<double, double>(0, 0)),
    _verbosity(-1)
{
  if ( getenv("limits_verbosity") != (char*)0 )
    _verbosity = atoi(getenv("limits_verbosity"));
  
  OBJ = this;
  normalize();
  RooArgList list(obs);
  for(size_t c=0; c < _data.size(); c++)
    {
      RooAbsArg* arg = list.at(c);
      _data[c] = (dynamic_cast<RooRealVar*>(arg))->getVal();
    }
}
#endif


Bayes::~Bayes() 
{
  // if _rfpoi is non-zero, this means we are using the RooFit
  // interface
  #ifdef __WITH_ROOFIT__
  if (_rfpoi) delete _pdf;
  #endif
  if (_interp) delete _interp;
}

double 
Bayes::prior(double poi)
{
  if (_prior)
    {
      return (*_prior)(poi);
    }
#ifdef __WITH_ROOFIT__
  else if (_rfprior)
    {
      _rfpoi->setVal(poi);
      return _rfprior->getVal();
    }
#endif
  return 1;
}

double 
Bayes::likelihood(double poi)
{
  return (*_pdf)(_data, poi);
}

double 
Bayes::normalize()
{
  // try to optimize support of likelihood x prior density

  int   nsteps = 2 * _nsteps;
  vector<double> p(nsteps+1);
  double step  = 0;
  
  for(int ii=0; ii < 2; ii++)
    {
      double pmax = 0.0;
      int    mode = 0;
      _poimax += step;      
      step = (_poimax - _poimin) / nsteps;
      
      for(int i=0; i < nsteps+1; i++)
  	{
  	  double xx = _poimin + i*step;
  	  p[i] = _likeprior(xx);
      
  	  if ( p[i] > pmax )
  	    {
  	      pmax = p[i];
  	      mode = i;
  	    }
  	}

      // find left edge of support
      double factor = 1.e-5;
      int jj = 0;
      for(int i=0; i < mode; i++)
  	{
  	  if (p[i] < factor * pmax)
  	    jj = i;
  	  else
  	    break;
  	}
      _poimin = jj * step;
  
      // find right edge of support  
      for(int i=mode; i < nsteps+1; i++)
  	{
  	  if (p[i] < factor * pmax)
  	    break;
  	  else
  	    jj = i;
  	}
      _poimax = 2 * jj * step;
    }
  assert( _poimax > _poimin );

  // now that wew have the support, calculate the unnormalized posterior
  // density at equal intervals;
  step = (_poimax - _poimin) / nsteps;
  for(int i=0; i < nsteps+1; i++)
    {
      double xx = _poimin + i*step;
      p[i] = _likeprior(xx);
    }

  
  // ROOT::Math::WrappedMemFunction<Bayes, double (Bayes::*)(double)> 
  //   fn(*this, &Bayes::_likeprior);
  // ROOT::Math::Integrator ifn(fn);
  // _normalization = ifn.Integral(_poimin, _poimax);
  // _normalize = false;

  // // compute posterior density over its support
  // step = (_poimax - _poimin)/ _nsteps;
  // for(int i=0; i < _nsteps+1; i++)
  //   {
  //     double poi = _poimin + i*step;
  //     p[i] = posterior(poi);
  //     //printf("%10.3f\t%10.3e\n", poi, p[i]);
  //   }

  // Compute cdf at several points

  step *= 2;
  double scale = step / 6;
  _x.resize(_nsteps+1);
  _y.resize(_nsteps+1);
  _x[0] = _poimin;
  _y[0] = 0;
  for(int i=1; i < nsteps+1; i+=2)
    {
      int j = (i + 1) / 2;
      _x[j] = _poimin + j * step;
      _y[j] = _y[j-1] + scale *  (p[i+1] + 4 * p[i] + p[i-1]);
      //printf("i[%d], j[%d], _x[%f], _y[%e]\n", i, j, _x[j], _y[j]);
    }
  assert(abs(_x.back() - _poimax) < 1.e-4 * _poimax);
  _poimax = _x.back();
  
  _normalization = _y.back();
  for(int i=0; i < _nsteps+1; i++) _y[i] /= _normalization;
  _normalize = false;

  if ( _interp != 0 )
    {
      delete _interp;
      _interp = 0;
    }
  if ( _interp == 0)
    _interp = new ROOT::Math::Interpolator(_x.size(),
					   ROOT::Math::
					   Interpolation::kLINEAR);
					   //Interpolation::kCSPLINE);
  _interp->SetData(_x, _y);

  return _normalization;
}

double 
Bayes::posterior(double poi)
{
  if ( _normalize ) normalize();
  return _likeprior(poi) / _normalization;
}

double 
Bayes::cdf(double poi)
{
  if ( poi < _poimin ) 
    return 0;
  else if ( poi > _poimax )
    return 1;

  if ( _normalize ) normalize();
  // Compute CDF
  ROOT::Math::WrappedMemFunction<Bayes, double (Bayes::*)(double)> 
    fn(*this, &Bayes::posterior);
  ROOT::Math::Integrator ifn(fn);
  return ifn.Integral(_poimin, poi);
}

double 
Bayes::percentile(double p)
{
  if ( p > 0 ) _cl = p; // Credibility level
  if ( _normalize ) normalize();

  // function whose root is to be found
  ROOT::Math::WrappedMemFunction<Bayes, double (Bayes::*)(double)> 
    fn(*this, &Bayes::_q);
  ROOT::Math::RootFinder rootfinder;
  
  rootfinder.SetFunction(fn, _poimin, _poimax);
  int status = rootfinder.Solve();
  if ( status != 1 )
    {
      cout << "*** Bayes *** RootFinder failed to find quantile"
           << endl;
      return -1;
    }
  return rootfinder.Root();
}

pair<double, double>
Bayes::MAP(double cl_)
{
  if ( _MAPdone ) return _result;
  if ( _normalize ) normalize();
  
  TMinuit minuit(1);
  minuit.SetPrintLevel(_verbosity);
  minuit.SetFCN(nlpFunc);

  double chi2 = TMath::NormQuantile((1.0+cl_)/2);
  chi2 = chi2*chi2/2;  
  minuit.SetErrorDef(chi2);

  int status=0;
  double guess = (_poimax+_poimin)/2;
  double stepsize = (_poimax-_poimin)/10000;
  minuit.mnparm(0, "poi", guess, stepsize, 
		_poimin, _poimax, status);
  double args[2] = {MAXITER, TOLERANCE};
  minuit.mnexcm("MIGRAD", args, 2, status);
  if ( status != 0 )
    {
      cout << "Bayes::MAP failed to find MAP" << endl;
      return pair<double, double>(0,0);
    }
  
  // get fit result
  minuit.GetParameter(0, _result.first, _result.second);
  
  _MAPdone = true;
  
  // // get MINOS errors
  // double upper, lower, gcc;
  // minuit.mnerrs(0, upper, lower, perror, gcc);
  // if ( (int)results.size() > 0 ) results[0] = poihat;
  // if ( (int)results.size() > 1 ) results[1] = upper;
  // if ( (int)results.size() > 2 ) results[2] = lower;
  // if ( (int)results.size() > 3 ) results[3] = poierr;
  // if ( (int)results.size() > 4 ) results[4] = gcc;
  return _result;
}

void
Bayes::setData(std::vector<double>& d)
{
  _data = d;
  normalize();
  _MAPdone = false;
}

double 
Bayes::zvalue(double poi)
{
  double N = likelihood(poi);
  double D = likelihood(0);
  double lnB10 = log(N/D);
  if ( lnB10 != lnB10 )
    {
      cout << "*** Bayes - this is Baaaad! lnB10 = " 
           << lnB10 << endl;
      exit(0);
    }  
  double abslnB10 = abs(lnB10);
  return lnB10 * sqrt(2*abslnB10) / abslnB10;
}

double
Bayes::estimate()
{
  MAP();
  return _result.first;
}

double
Bayes::uncertainty()
{
  MAP();
  return _result.second;
}

double 
Bayes::_likeprior(double poi)
{
  if ( poi != poi )
    {
      cout << "*** Bayes - this is Baaaad! poi = " 
           << poi << endl;
      exit(0);
    }
  return likelihood(poi) * prior(poi);
}

double 
Bayes::_q(double poi)
{
  if ( poi < _poimin ) 
    return 0;
  else if ( poi > _poimax )
    return 1;
  else
     return _interp->Eval(poi) - _cl;
}

