//--------------------------------------------------------------
// File: Bayes.cc
// Description: Implement standalone Bayes calculator
// 
// Created: 11 Jan 2011 HBP
// Updated: 13 Mar 2011 HBP - fix normalization of post.
//          06 Mar 2014 HBP
//--------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
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
  const int MAXITER=1000;
  const double TOLERANCE=1.e-4;
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
	     double cl,
	     PriorFunction* prior_)
  : _pdf(&model),
    _data(d),
    _poimin(poimin),
    _poimax(poimax),
    _cl(cl),
    _prior(prior_),
    _normalize(true),
    _interp(0),
    _nsteps(200),
    _x(vector<double>()),
    _y(vector<double>())
{
  OBJ = this;
}


Bayes::~Bayes() 
{
  if (_interp) delete _interp;
}

double 
Bayes::prior(double poi)
{
  if (_prior)
    {
      return (*_prior)(poi);
    }
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
  ROOT::Math::WrappedMemFunction<Bayes, double (Bayes::*)(double)> 
    fn(*this, &Bayes::_likeprior);
  ROOT::Math::Integrator ifn(fn);
  _normalization = ifn.Integral(_poimin, _poimax);
  _normalize = false;
  
  // Compute cdf at several points
  vector<double> p(_nsteps+1);
  double step = (_poimax - _poimin)/_nsteps;
  for(int i=0; i < _nsteps+1; i++)
    {
      double xx = _poimin + i*step;
      p[i] = posterior(xx);
    }

  vector<double> c1(_nsteps, 0);
  for(int i=1; i < _nsteps+1; i++)
    c1[i] = c1[i-1] + 0.5*step*(p[i]+p[i-1]);
  
  _x.clear();
  _y.clear();
  double c2 = 0;
  _x.push_back(_poimin);
  _y.push_back(c2);
  for(int i=2; i < _nsteps+1; i+=2)
    {
      double p1 = p[i-2];
      double p2 = p[i];
      c2 += step*(p2+p1);
      
      _x.push_back(_poimin + i*step);          
      double z = (4*c1[i]-c2)/3;
      if ( z > 1 ) z = 1;
      _y.push_back(z);
    }

  if ( _interp == 0)
    _interp = new ROOT::Math::Interpolator(_x.size(),
					   ROOT::Math::
					   Interpolation::kCSPLINE);
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
Bayes::quantile(double p)
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

void
Bayes::MAP(double CL, vector<double>& results)
{
  if ( _normalize ) normalize();

  TMinuit minuit(1);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(nlpFunc);

  double chi2 = TMath::NormQuantile((1.0+CL)/2);
  chi2 = chi2*chi2/2;
  minuit.SetErrorDef(chi2);

  int status=0;
  double guess = (_poimax+_poimin)/2;
  double stepsize = (_poimax-_poimin)/100;
  minuit.mnparm(0, "poi", guess, stepsize, 
		_poimin, _poimax, status);
  double args[2] = {MAXITER, TOLERANCE};
  minuit.mnexcm("MIGRAD", args, 2, status);
  if ( status != 0 )
    {
      cout << "Bayes::MAP failed to find MAP" << endl;
      return;
    }
  
  // get fit result
  double poihat = 0.0;
  double poierr = 0.0;
  minuit.GetParameter(0, poihat, poierr);

  // get MINOS errors
  double upper, lower, perror, gcc;
  minuit.mnerrs(0, upper, lower, perror, gcc);
  if ( (int)results.size() > 0 ) results[0] = poihat;
  if ( (int)results.size() > 1 ) results[1] = upper;
  if ( (int)results.size() > 2 ) results[2] = lower;
  if ( (int)results.size() > 3 ) results[3] = perror;
  if ( (int)results.size() > 4 ) results[4] = gcc;
  return;
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

