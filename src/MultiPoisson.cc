//--------------------------------------------------------------
// File: MultiPoisson.cc
// Description: Implement the multi-Poisson model 
// 
// Created: June 11, 2014
//--------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <cmath>
#include "TMath.h"
#include "MultiPoisson.h"

using namespace std;

MultiPoisson::MultiPoisson()
  : PDFunction(),
    _N(vector<double>()),
    _Ngen(vector<double>()),
    _L(vector<double>()),
    _eff(vector<vector<double> >()),
    _bkg(vector<vector<double> >()),
    _meanefl(vector<double>()),
    _meanbkg(vector<double>()),    
    _random(TRandom3()),
    _nbins(0),
    _index(-1),
    _profile(false)
{}

MultiPoisson::MultiPoisson(vector<double>& N)
  : PDFunction(),
    _N(N),
    _Ngen(N),
    _L(vector<double>()),
    _eff(vector<vector<double> >()),
    _bkg(vector<vector<double> >()),
    _meanefl(vector<double>(N.size(),0)),
    _meanbkg(vector<double>(N.size(),0)),      
    _random(TRandom3()),
    _nbins((int)N.size()),
    _index(-1),
    _profile(false)
{}

MultiPoisson::~MultiPoisson() 
{}

void MultiPoisson::add(double L, vector<double>& eff, vector<double>& bkg)
{
  _L.push_back(L);
  _eff.push_back(eff);
  _bkg.push_back(bkg);
}

void MultiPoisson::update(int ii, vector<double>& eff)
{
  if ( ii < 0 ) return;
  if ( ii > (int)(_eff.size()-1) ) return;
  copy(eff.begin(), eff.end(), _eff[ii].begin());
}

void MultiPoisson::computeMeans()
{
  int N = _L.size();
  for(int ibin=0; ibin < _nbins; ++ibin)
    {
      _meanefl[ibin] = 0.0;
      _meanbkg[ibin] = 0.0;
      for(int ii=0; ii < N; ++ii)
	{
	  _meanefl[ibin] += _eff[ii][ibin]*_L[ii];
	  _meanbkg[ibin] += _bkg[ii][ibin];
	}
      _meanefl[ibin] /= N;
      _meanbkg[ibin] /= N;
    }  
}

void MultiPoisson::set(int ii)
{
  if ( ii < 0 ) return;
  if ( ii > (int)(_eff.size()-1) ) return;
  _index = ii;
}

void MultiPoisson::reset()
{
  _index = -1;
}

vector<double>
MultiPoisson::get(int ii)
{
  if ( ii < 0 ) return vector<double>();
  if ( ii > (int)(_eff.size()-1) ) return vector<double>();
  vector<double> d(1+2*_eff[ii].size());
  d[0] = _L[ii];
  double* a = &d[0];
  copy(_eff[ii].begin(), _eff[ii].end(), a+1);
  copy(_bkg[ii].begin(), _bkg[ii].end(), a+1+_eff[ii].size());
  return d;
}

vector<double>&  
MultiPoisson::generate(double sigma)
{
  if(_nbins == 0)
    {
      cout << "MultiPoisson::generate: nbins = 0" << endl;
      exit(0);
    }

  int nconstants = _L.size();
  int icon = _random.Integer(nconstants-1);      
  for(int ibin=0; ibin < _nbins; ++ibin)
    {
      double mean = _eff[icon][ibin] * _L[icon] * sigma + _bkg[icon][ibin];
      _Ngen[ibin] = _random.Poisson(mean);
    }
  return _Ngen;
}

vector<double>&  
MultiPoisson::generateUsingMeans(double sigma)
{
  if(_nbins == 0)
    {
      cout << "MultiPoisson::generate: nbins = 0" << endl;
      exit(0);
    }
  for(int ibin=0; ibin < _nbins; ++ibin)
    {
      double mean = _meanefl[ibin] * sigma + _meanbkg[ibin];
      _Ngen[ibin] = _random.Poisson(mean);
    }
  return _Ngen;
}

double 
MultiPoisson::operator() (std::vector<double>& N, double sigma)
{
  int first = 0;
  int last  = _L.size()-1;
  if ( _index >= 0 )
    {
      first = _index;
      last  = _index;
    }
  int nconstants = 1 + last - first;

  double likelihood = 0.0;
  if ( _profile )
    {
      // do something!
    }
  else
    {
      for(int icon=first; icon <= last; ++icon)
	{
	  double p = 1.0;
	  for(int ibin=0; ibin < _nbins; ++ibin)
	    {
	      double mean = _eff[icon][ibin]*_L[icon]*sigma + _bkg[icon][ibin];
	      p *= TMath::Poisson(N[ibin], mean);
	    }
	  likelihood += p;
	}
      likelihood /= nconstants;
    }
  return likelihood;
}

void 
MultiPoisson::setSeed(int seed) { _random.SetSeed(seed); }

double 
MultiPoisson::operator() (double sigma)
{
  return (*this)(_N, sigma);
}
