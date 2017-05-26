//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: MultiPoissonGammaModel.cc
// Description: Implement the Poisson-Gamma model 
//              marginalized over the nuisance parameters.
// 
// Created: 11-Jun-2010
// Updated: 04-Jul-2015 HBP Renamed MultiPoissonGammaModel.cc 
//
//--------------------------------------------------------------
#include <iostream>
#include <cmath>
#include "TError.h"
#include "MultiPoissonGammaModel.h"

using namespace std;
//--------------------------------------------------------------
MultiPoissonGammaModel::MultiPoissonGammaModel()
  : PDFunction(),
    _data(vector<double>()),
    _x(vector<double>()),
    _a(vector<double>()),
    _y(vector<double>()),
    _b(vector<double>()),
    _maxcount(100000),
    _gslRan(new ROOT::Math::Random<ROOT::Math::GSLRngMT>())
{}

MultiPoissonGammaModel::MultiPoissonGammaModel(vector<double>& data,
					       vector<double>& x,
					       double a,
					       vector<double>& y,
					       double b, 
					       int maxcount)
  : PDFunction(),
    _data(data),
    _x(x),
    _a(vector<double>(x.size(), a)),
    _y(y),
    _b(vector<double>(x.size(), b)),
    _maxcount(maxcount),
    _gslRan(new ROOT::Math::Random<ROOT::Math::GSLRngMT>())
{
}

MultiPoissonGammaModel::MultiPoissonGammaModel(vector<double>& data,
					       vector<double>& x, 
					       vector<double>& a,
					       vector<double>& y, 
					       vector<double>& b,
					       int maxcount)
  : PDFunction(),
    _data(data),
    _x(x),
    _a(a),
    _y(y),
    _b(b),
    _maxcount(maxcount),
    _gslRan(new ROOT::Math::Random<ROOT::Math::GSLRngMT>())
{
  if(_x.size() != _b.size() ||
     _x.size() != _a.size() ||
     _a.size() != _b.size())
    {
      Error("MultiPoissonGammaModel",
	    "input vectors have different sizes.");
      exit(0);
    }
}

MultiPoissonGammaModel::MultiPoissonGammaModel(double data,
					       double x,
					       double a,
					       double y,
					       double b, 
					       int maxcount)
  : PDFunction(),
    _data(vector<double>(1, data)),
    _x(vector<double>(1, x)),
    _a(vector<double>(1, a)),
    _y(vector<double>(1, y)),
    _b(vector<double>(1, b)),
    _maxcount(maxcount),
    _gslRan(new ROOT::Math::Random<ROOT::Math::GSLRngMT>())
{
}

MultiPoissonGammaModel::~MultiPoissonGammaModel() 
{}

vector<double>&  
MultiPoissonGammaModel::generate(double sigma)
{
  if((int)_x.size() == 0)
    {
      Error("MultiPoissonGammaModel",
	    "required input vectors not supplied by user.");
      exit(0);
    }
      
  for(size_t ibin=0; ibin < _x.size(); ++ibin)
    {
      double epsilon = _gslRan->Gamma(_x[ibin]+0.5, 1.0/_a[ibin]);
      //cout << "epsilon " << epsilon << endl;
            
      double mu      = _gslRan->Gamma(_y[ibin]+0.5, 1.0/_b[ibin]);
      //cout << "mu      " << mu << endl;
      
      double mean    = epsilon * sigma + mu;
      //cout << "mean    " << mean << endl;
      
      _data[ibin]    = _gslRan->Poisson(mean);
      //cout << "data  = " << _data[ibin] << endl;
    }
  return _data;
}

double 
MultiPoissonGammaModel::operator() (std::vector<double>& data, double sigma)
{
  if(data.size() != _x.size())
    {
      Error("MultiPoissonGammaModel",
	    "input vector size != %d bins", (int)_x.size());
      exit(0);
    }
  long double C1[_maxcount+1];
  long double C2[_maxcount+1];
    
  // loop over bins
  long double prob = 1.0;    
  for(size_t ibin=0; ibin < _x.size(); ++ibin)
    {
      double nn = data[ibin];    // observed count	  
      if ( nn > _maxcount )
	{
          Error("MultiPoissonGammaModel",
                "bin %d has a count, %d, that "
		"is greater than maxcount, %d.",
                (int)ibin, (int)nn, _maxcount);
          exit(0);
	}
      
      double p1 = sigma / _a[ibin];
      double p2 = 1.0 / _b[ibin];
      double A1 = _x[ibin]-0.5;  // signal count
      double A2 = _y[ibin]-0.5;  // background count

      // compute coefficients
      C1[0] = pow(1+p1, -(A1+1));
      C2[0] = pow(1+p2, -(A2+1));
      double dk;	  
      if ( nn > 0 )
	{
          for(int ik=1; ik <= (int)nn; ++ik)
	    {
              dk = (double)ik;
              C1[ik] = C1[ik-1] * (p1/(1+p1)) * (A1+dk)/dk;
              C2[ik] = C2[ik-1] * (p2/(1+p2)) * (A2+dk)/dk;
	    }
	}
	  
      // compute p(nn|sigma)
      long double sum = 0.0;  
      for (int ik=0; ik <= (int)nn; ++ik)
	{
          sum += C1[ik] * C2[(int)nn-ik];
	}	  
      prob *= sum; // product of likelihood over bins	  
    } // loop over bins
  return (double)prob;
}

double 
MultiPoissonGammaModel::operator() (double sigma)
{
  return (*this)(_data, sigma);
}
