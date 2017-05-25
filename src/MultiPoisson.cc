//--------------------------------------------------------------
// File: MultiPoisson.cc
// Description: Implement the multi-Poisson model 
// 
// Created: 11-Jun-2014 HBP
//          25-May-2017 HBP simplify config file format
//                      warning: not backwards compatible!
//--------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include "TMath.h"
#include "MultiPoisson.h"
#include "TError.h"

using namespace std;

namespace {
  ///
  std::string strip(std::string line)
  {
    int l = line.size();
    if ( l == 0 ) return std::string("");
    int n = 0;
    while (((line[n] == 0)    ||
	    (line[n] == ' ' ) ||
	    (line[n] == '\n') ||
	    (line[n] == '\t')) && n < l) n++;
    
    int m = l-1;
    while (((line[m] == 0)    ||
	    (line[m] == ' ')  ||
	    (line[m] == '\n') ||
	    (line[m] == '\t')) && m > 0) m--;
    return line.substr(n,m-n+1);
  }
};


MultiPoisson::MultiPoisson()
  : PDFunction(),
    _N(vector<double>()),
    _Ngen(vector<double>()),
    _S(vector<vector<double> >()),
    _B(vector<vector<double> >()),
    _meanS(vector<double>()),
    _meanB(vector<double>()),    
    _random(TRandom3()),
    _nbins(0),
    _index(-1),
    _profile(false)
{}


MultiPoisson::MultiPoisson(string filename)
  : PDFunction(),
    _N(vector<double>()),
    _Ngen(vector<double>()),
    _S(vector<vector<double> >()),
    _B(vector<vector<double> >()),
    _meanS(vector<double>()),
    _meanB(vector<double>()),    
    _random(TRandom3()),
    _nbins(0),
    _index(-1),
    _profile(false) 
{
  // open input text file with format
  // number of bins  ... 
  // count1 count2 ...
  // number of points
  // S1   S2   ...
  // B1   B2   ...
  
  ifstream inp(filename.c_str());
  if ( ! inp.good() )
    {
      Error("MultiPoisson", "unable to open file %s", filename.c_str());
      exit(0);
    }
  
  // read lines
  vector<string> records;
  string line;
  while(getline(inp, line))
    {
      line = strip(line);
      if (line == "") continue;
      if (line.substr(0,1) == "#") continue;
      records.push_back(line);
    }

  if (records.size() == 0)
    {
      Error("MultiPoisson", "zero uncommented records in file %s",
	    filename.c_str());
      exit(0);
    }
  
  // get number of bins
  _nbins = 0;
  int lineno = 0;
  stringstream nin(records[lineno]);
  try
    {
      nin >> _nbins;
    }
  catch (...)
    {
      Error("MultiPoisson",
	    "unable to get number of bins; check file format");
      exit(0);
    }

  // get observed counts 
  
  vector<double> S;
  vector<double> B;
  lineno++;
  stringstream din(records[lineno]);

  int c = 0;
  double x;
  for(int i=0; i < _nbins; i++)
    {
      try
	{
	  c++;
	  din >> x;
	}
      catch (...)
	{
	  Error("MultiPoisson", "problem accessing observed count at bin %d",
		c);
	  exit(0);	  
	}
      _N.push_back(x);
      _meanS.push_back(0);
      _meanB.push_back(0);
      S.push_back(0);
      B.push_back(0);
    }
  
  // get number of samples
  int samplesize = 0;
  lineno++;
  stringstream hin(records[lineno]);
  try
    {
      hin >> samplesize;
    }
  catch (...)
    {
      Error("MultiPoisson",
	    "unable to get sample size; check file format");
      exit(0);
    }

  // loop over sampled points  
  for(int ii=0; ii < samplesize; ii++)
    {
      // get signals or effective luminosities
      lineno++;
      stringstream sin(records[lineno]);
      try
	{
	  for(int i=0; i < _nbins; i++) sin >> S[i];
	}
      catch (...)
	{
	  Error("MultiPoisson",
		"problem decoding signal counts:\n%s",
		records[lineno].c_str());
	  exit(0);
	}

      // get backgrounds
      lineno++;
      stringstream bin(records[lineno]);
      try
	{
	  for(int i=0; i < _nbins; i++) bin >> B[i];
	}
      catch(...)
	{
	  Error("MultiPoisson",
		"problem decoding background counts:\n%s",
		records[lineno].c_str());
	  exit(0);
	}

      add(S, B);
    }
  computeMeans();
}

MultiPoisson::MultiPoisson(vector<double>& N)
  : PDFunction(),
    _N(N),
    _Ngen(N),
    _S(vector<vector<double> >()),
    _B(vector<vector<double> >()),
    _meanS(vector<double>(N.size(),0)),
    _meanB(vector<double>(N.size(),0)),      
    _random(TRandom3()),
    _nbins((int)N.size()),
    _index(-1),
    _profile(false)
{}

MultiPoisson::~MultiPoisson() 
{}

 void MultiPoisson::add(vector<double>& S, vector<double>& B)
{
  _S.push_back(S);
  _B.push_back(B);
}

void MultiPoisson::update(int ii, vector<double>& S)
{
  if ( ii < 0 ) return;
  if ( ii > (int)(_S.size()-1) ) return;
  copy(S.begin(), S.end(), _S[ii].begin());
}

void MultiPoisson::computeMeans()
{
  int M = _S.size();
  for(int ibin=0; ibin < _nbins; ++ibin)
    {
      _meanS[ibin] = 0.0;
      _meanB[ibin] = 0.0;
      for(int ii=0; ii < M; ++ii)
	{
	  _meanS[ibin] += _S[ii][ibin];
	  _meanB[ibin] += _B[ii][ibin];
	}
      _meanS[ibin] /= M;
      _meanB[ibin] /= M;
    }  
}

void MultiPoisson::set(int ii)
{
  if ( ii < 0 ) return;
  if ( ii > (int)(_S.size()-1) ) return;
  _index = ii;
}

void MultiPoisson::reset()
{
  _index = -1;
}

// vector<pair<double, double> >
// MultiPoisson::get(int ii)
// {
//   vector<pair<double, double> > p;
//   if ( ii < 0 ) return p;
//   if ( ii > (int)(_S.size()-1) ) return p;

//   for(size_t c=0; c < _S.size(); c++)
//     p.push_back(pair<double, double>(_S[c], _B[c]));
//   return p;
// }

vector<double>&  
MultiPoisson::generate(double mu)
{
  if(_nbins <= 0)
    {
      Error("MultiPoisson", "nbins = 0, can't generate!");
      exit(0);
    }
  int nconstants = _S.size();
  int icon = _random.Integer(nconstants-1);      
  for(int ibin=0; ibin < _nbins; ++ibin)
    {
      double mean = mu * _S[icon][ibin] + _B[icon][ibin];
      _Ngen[ibin] = _random.Poisson(mean);
    }
  return _Ngen;
}

double 
MultiPoisson::operator() (std::vector<double>& N, double mu)
{
  int first = 0;
  int last  = _S.size()-1;
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
	      double mean = mu * _S[icon][ibin] + _B[icon][ibin];
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
MultiPoisson::operator() (double mu)
{
  return (*this)(_N, mu);
}
