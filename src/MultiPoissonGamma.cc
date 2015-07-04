//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: MultiPoissonGamma.cc
// Description: Implements the Poisson-Gamma model 
//              that is marginalized over the nuisance parameters.
// 
// Created: 11-Jun-2010
// Updated: 04-Jul-2015 HBP Renamed MultiPoissonGamma.cc 
//
//--------------------------------------------------------------
#include <iostream>
#include <cmath>
#include "TError.h"
#include "MultiPoissonGamma.h"

using namespace std;
//--------------------------------------------------------------
MultiPoissonGamma::MultiPoissonGamma()
  : PDFunction(),
    _yy(vector<double>()),
    _xx(vector<double>()),
    _nbins(0),
    _b(vector<double>()),
    _a(vector<double>()),
    _maxcount(100000),
    _xdata(vector<double>()),
    _gslRan(ROOT::Math::Random<ROOT::Math::GSLRngMT>())
{}

MultiPoissonGamma::MultiPoissonGamma(vector<double>& yy, 
				     vector<double>& xx,
				     double b, 
				     double a,
				     int maxcount)
  : PDFunction(),
    _yy(yy),
    _xx(xx),
    _nbins((int)yy.size()),
    _b(vector<double>(_nbins, b)),
    _a(vector<double>(_nbins, a)),
    _maxcount(maxcount),
    _xdata(vector<double>(_nbins, 0)),
    _gslRan(ROOT::Math::Random<ROOT::Math::GSLRngMT>())
{
  if((int)_xx.size() != _nbins)
    {
      Error("MultiPoissonGamma",
	    "input vectors have different sizes.");
      exit(0);
    }
}

MultiPoissonGamma::MultiPoissonGamma(vector<double>& yy, 
				     vector<double>& xx,
				     vector<double>& b, 
				     vector<double>& a,
				     int maxcount)
  : PDFunction(),
    _yy(yy),
    _xx(xx),
    _nbins((int)yy.size()),
    _b(b),
    _a(a),
    _maxcount(maxcount),
    _xdata(vector<double>(_nbins, 0)),
    _gslRan(ROOT::Math::Random<ROOT::Math::GSLRngMT>())
{
  if((int)_xx.size() != _nbins ||
     (int)_b.size()  != _nbins ||
     (int)_a.size()  != _nbins)
    {
      Error("MultiPoissonGamma",
	    "input vectors have different sizes.");
      exit(0);
    }
}

MultiPoissonGamma::~MultiPoissonGamma() 
{}

vector<double>&  
MultiPoissonGamma::generate(double sigma)
{
  if(_nbins == 0)
    {
      Error("MultiPoissonGamma",
	    "required input vectors not supplied by user.");
      exit(0);
    }
      
  for(int ibin=0; ibin < _nbins; ++ibin)
    {
      double epsilon = _gslRan.Gamma(_xx[ibin]+0.5, 1.0/_a[ibin]);
      double b       = _gslRan.Gamma(_yy[ibin]+0.5, 1.0/_b[ibin]);
      double mean    = epsilon * sigma + b;
      _xdata[ibin]   = _gslRan.Poisson(mean);
    }
  return _xdata;
}

double 
MultiPoissonGamma::operator() (std::vector<double>& xdata, double sigma)
{
  if((int)xdata.size() != _nbins)
    {
      Error("MultiPoissonGamma",
	    "input vector size != %d bins", _nbins);
      exit(0);
    }
  
  long double C1[_maxcount+1];
  long double C2[_maxcount+1];
    
  // loop over bins
  long double prob = 1.0;    
  for(int ibin=0; ibin < _nbins; ++ibin)
    {
      double nn = xdata[ibin];    // observed count	  
      if ( nn > _maxcount )
	{
          Error("MultiPoissonGamma",
                "bin %d has a count: %10.3e that is greater than maxcount: %d.", 
                ibin, nn, _maxcount);
          exit(0);
	}
      
      double p1 = sigma / _a[ibin];
      double p2 = 1.0 / _b[ibin];
      double A1 = _xx[ibin]-0.5;  // signal count
      double A2 = _yy[ibin]-0.5;  // background count

      // compute coefficients
      C1[0] = pow(1+p1, -(A1+1));
      C2[0] = pow(1+p2, -(A2+1));
      double dk;	  
      if ( nn > 0 )
	{
          for(int ik=1; ik <= (int)nn; ++ik)
	    {
              dk = (double)ik;
              C1[ik] = C1[ik-1]*(p1/(1+p1))*(A1+dk)/dk;
              C2[ik] = C2[ik-1]*(p2/(1+p2))*(A2+dk)/dk;
	    }
	}
	  
      // compute p(nn|sigma)
      long double sum = 0.0;  
      for (int ik=0; ik <= (int)nn; ++ik)
	{
          sum += C1[ik]*C2[(int)nn-ik];
	}	  
      prob *= sum; // product of likelihood over bins	  
    } // loop over bins
  return (double)prob;
}

