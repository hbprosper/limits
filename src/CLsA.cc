// ---------------------------------------------------------------------------
// File: CLsA.cc
// Description: compute CLsA limits using statistic
//              Q = -2*ln L(poi)/L(poi_hat), where poi_hat >= 0
// 
// Created: 07 Sep 2011 HBP
// Updated: 05 Mar 2014 HBP 
//          11 Aug 2014 HBP - fix bracketing bug
//          13 Aug 2014 HBP - cache distributions for reuse
//          14 Aug 2014 HBP - use TH1 to smooth and interpolate
//          03 Jun 2015 HBP - divorce CLsA from CLs
// ---------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include "TH1D.h"
#include "TMinuit.h"
#include "TMath.h"
#include "Math/WrappedFunction.h"
#include "Math/RootFinder.h"
#include "CLsA.h"

ClassImp(CLsA);

using namespace std;
// ---------------------------------------------------------------------------
// function to be minimized
namespace {
  const int MAXITER=1000;
  const double TOLERANCE=1.e-3;
  CLsA* OBJ=0;
  void nllFunc(int&    /*npar*/, 
	       double* /*grad*/, 
	       double& fval,
	       double* xval,
	       int     /*iflag*/)
  {
    double poi = xval[0];
    fval = OBJ->nll(poi);
  }
};

CLsA::CLsA(PDFunction& model,
	   vector<double>& data,
	   double poimin,   // minimum value of parameter of interest
	   double poimax,   // maximum value of parameter of interest
	   double CL)       // alpha = 1 - CL
  : _model(&model),
    _data(data),
    _poimin(poimin),
    _poimax(poimax),
    _alpha(1-CL),
    _CLsbA(0),
    _CLbA(0),
    _Qobs(0),
    _Ddata(0)
{
  OBJ = this;
  // find best fit using specified data
  fit();
  // cache denominator computed using data
  _Ddata = (*_model)(_data, _poihat);
}

CLsA::~CLsA() {}

void CLsA::setData(std::vector<double>& d) 
{ 
  _data = d;
  fit();
  // cache denominator computed using data
  _Ddata = (*_model)(_data, _poihat); 
}

double CLsA::fit(double guess)
{
  _poihat = 0.0;
  _poierr = 0.0;
  
  TMinuit minuit(1);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(nllFunc);
  minuit.SetErrorDef(0.5); // for log-likelihood fit

  int status=0;
  if ( guess < 0 ) guess = (_poimax+_poimin)/2;
  double stepsize = (_poimax-_poimin)/100;
  minuit.mnparm(0, "poi", guess, stepsize, 
		_poimin, _poimax, status);
  double args[2] = {MAXITER, TOLERANCE};
  minuit.mnexcm("MIGRAD", args, 2, status);
  if ( status != 0 )
    {
      cout << "CLsA::fit failed to find MLE" << endl;
      cout << "CLsA::fit status  = " << status    << endl;
      cout << "CLsA::fit _poimin = " << _poimin  << endl;
      cout << "CLsA::fit _poimax = " << _poimax  << endl;
      cout << "CLsA::fit stepsize= " << stepsize << endl;
      cout << "CLsA::fit guess   = " << guess    << endl;
      _poihat = 0;
      cout << "CLsA::fit - set _poihat = " << _poihat << endl;
      return _poihat;
    }

  // get fit result
  minuit.GetParameter(0, _poihat, _poierr);

  return _poihat;
}
 
double CLsA::nll(double poi)
{
  return -log((*_model)(_data, poi));
}

double CLsA::Qstatistic(double poi)
{
  // Compute statistic
  double N = (*_model)(_data,  poi);
  if ( N != N )
    {
      cout << "CLsA::Qstatistic - N is Nan for poi = " << poi 
	   << ", setting Q = 0.0" << endl;
      return 0.0;
    }

  double D = _Ddata;
  if ( D != D )
    {
      cout << "CLsA::Qstatistic - D is Nan for poi = " << poi 
	   << ", setting Q = 0.0" << endl;
      return 0.0;
    }

  try
    {
      double R = N / D;
      if ( R != R )
	{
	  cout << "CLsA::Qstatistic - R = N/D is Nan for poi = " << poi 
	       << ", setting Q = 0.0" << endl;
	  return 0.0;	  
	}

      // should always be positive, but because of rounding errors
      // Q could become very slightly negative
      double Q = -2*log(R);
      if ( Q > 0 )
	return Q;
      else
	return 0.0;
    }
  catch (...)
    {
      return 0;
    }
}

double CLsA::operator()(double poi)
{
  // compute observed value of statistic for given
  // parameter of interest (poi)
  _Qobs = Qstatistic(poi);
  if ( _Qobs != _Qobs )
    {
      cout << "CLsA::operator() - _Qobs is Nan at poi = " << poi << endl;
      cout << "CLsA::operator() - setting _Qobs = 0" << endl;
      _Qobs = 0.0;
     }

  // compute standard deviation at best fit value
  double h  = 1.e-3;
  double f0 = nll(_poihat);
  double f1 = nll(_poihat+h);
  double f2 = nll(_poihat+2*h);
  double d2 = (f0 - 2*f1 + f2)/(h*h); // 2nd derivative
  double stdv = 1.0/sqrt(fabs(d2));
  double z = sqrt(_Qobs);
  double y = z;
      
  // compute CL(s+b)
  _CLsbA = 1 - TMath::Freq(z);

  // compute CL(b)
  _CLbA  = TMath::Freq(poi/stdv-y);

  if ( _CLbA > 0 )
    return _CLsbA / _CLbA;
  else
    return -1;
}

double CLsA::limit(double CL)
{
  if ( CL > 0 ) _alpha = 1-CL;

  // function whose root is to be found
  ROOT::Math::WrappedMemFunction<CLsA, double (CLsA::*)(double)> 
    fn(*this, &CLsA::_f);
  ROOT::Math::RootFinder rootfinder;
  rootfinder.SetFunction(fn, _poimin, _poimax);
  int status = rootfinder.Solve();
  if ( status != 1 )
    {
      cout << "** CLsA::limit - RootFinder failed"
	   << endl;
    }
  return rootfinder.Root();
}

double CLsA::_f(double poi)
{
  return (*this)(poi) - _alpha;
}
