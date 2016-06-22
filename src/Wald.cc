// ---------------------------------------------------------------------------
// File: Wald.cc
// Description: compute limits using the statistic
//
//              q(poi) = -2*ln L(poi)/L(poi_hat)
//
//              and assuming that the Wald approximation
//
//              q(poi) = (poi-poi_hat)^2 / var(poi_hat), poi_hat < poi
//                    or 0 for poi_hat > poi.
//
//  See "Asymptotic formulae for likelihood-based tests of new physics",
//      G. Cowan, K. Cranmer, E. Gross, and O. Vitells, arXiv:1007.1727v3
// 
// Created: 04-Jun-2015 Harrison B. Prosper
// ---------------------------------------------------------------------------
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "TMinuit.h"
#include "TMath.h"
#include "Math/WrappedFunction.h"
#include "Math/RootFinder.h"
#include "Wald.h"

ClassImp(Wald);

using namespace std;
// ---------------------------------------------------------------------------
// function to be minimized
namespace {
  const int MAXITER=1000;
  const double TOLERANCE=1.e-3;
  Wald* OBJ=0;
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

Wald::Wald(PDFunction& model,
	   vector<double>& data,
	   double poimin,   // minimum value of parameter of interest
	   double poimax,   // maximum value of parameter of interest
	   double CL)       // alpha = 1 - CL
  : _model(&model),
    _data(data),
    _poimin(poimin),
    _poimax(poimax),
    _alpha(1-CL)
{
  OBJ = this;
  // find best fit value of parameter of interest
  fit();
}

Wald::~Wald() {}

void Wald::setData(std::vector<double>& d) 
{ 
  _data = d;
  fit();
}

double Wald::fit(double guess)
{
  _poihat = 0.0;
  _poierr = 0.0;
  
  TMinuit minuit(1);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(nllFunc);
  minuit.SetErrorDef(0.5); // for log-likelihood fit

  int status=0;
  if ( guess < 0 ) guess = (_poimax + _poimin)/2;
  double stepsize = (_poimax - _poimin)/100;
  minuit.mnparm(0, "poi", guess, stepsize, 
		_poimin, _poimax, status);
  double args[2] = {MAXITER, TOLERANCE};
  minuit.mnexcm("MIGRAD", args, 2, status);
  if ( status != 0 )
    {
      _poihat = 0;
      cout << "Wald::fit failed to find MLE" << endl;
      cout << "Wald::fit status  = " << status    << endl;
      cout << "Wald::fit _poimin = " << _poimin  << endl;
      cout << "Wald::fit _poimax = " << _poimax  << endl;
      cout << "Wald::fit stepsize= " << stepsize << endl;
      cout << "Wald::fit guess   = " << guess    << endl;
      cout << "Wald::fit - set _poihat = " << _poihat << endl;
      return _poihat;
    }
  // get fit result
  minuit.GetParameter(0, _poihat, _poierr);
  return _poihat;
}
 
double Wald::estimate()
{
  return _poihat;
}

double Wald::uncertainty()
{
  return _poierr;
}

double Wald::nll(double poi)
{
  return -log((*_model)(_data, poi));
}

double Wald::operator()(double poi)
{
  // compute observed value of statistic
  double qobs = 0;
  qobs = 2*(nll(poi) - nll(_poihat));
  if ( qobs != qobs )
    {
      cout << "Wald::operator() - qobs is Nan at poi = " << poi << endl;
      cout << "Wald::operator() - setting qobs = 0" << endl;
      qobs = 0.0;
     }
  // compute p(poi)
  double qobsabs = abs(qobs);
  double sign = 1.0;
  if ( qobs != 0 ) sign = qobs/qobsabs;
  double Z = sign*sqrt(abs(qobs));
  return 1 - TMath::Freq(Z);
}

double Wald::zvalue(double poi)
{
  // compute observed value of statistic
  double qobs = 0;
  qobs = 2*(nll(poi) - nll(_poihat));
  if ( qobs != qobs )
    {
      cout << "Wald::zvalue(poi) - qobs is Nan at poi = " << poi << endl;
      cout << "Wald::zvalue(poi) - setting qobs = 0" << endl;
      qobs = 0.0;
     }
  if ( qobs > 0 )
    return sqrt(qobs);
  else
    {
      cout << "Wald::zvalue(poi) - qobs = "
	   << qobs << " < 0  of "
	   << poi << endl;
      cout << "Wald::zvalue(poi) - setting Z = sign*sqrt(abs(qobs))"
	   << endl;      
      return -sqrt(abs(qobs));
    }
}


double Wald::quantile(double CL)
{
  if ( CL > 0 ) _alpha = 1-CL;

  // function whose root is to be found
  ROOT::Math::WrappedMemFunction<Wald, double (Wald::*)(double)> 
    fn(*this, &Wald::_f);
  ROOT::Math::RootFinder rootfinder;
  double poimin = _poimin;
  double poimax = _poimax;
  double alphapoi = (*this)(_poihat);
  
  if ( alphapoi > _alpha )
    poimin = _poihat;
  else
    {
      _alpha = CL;
      poimax = _poihat;
    }
  
  rootfinder.SetFunction(fn, poimin, poimax);
  int status = rootfinder.Solve();
  if ( status != 1 )
    {
      cout << "** Wald::limit - RootFinder failed"
	   << endl;
    }
  return rootfinder.Root();
}

double Wald::_f(double poi)
{
  return (*this)(poi) - _alpha;
}
