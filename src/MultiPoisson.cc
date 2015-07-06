//--------------------------------------------------------------
// File: MultiPoisson.cc
// Description: Implement the multi-Poisson model 
// 
// Created: June 11, 2014
//--------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "TMath.h"
#include "MultiPoisson.h"

using namespace std;

MultiPoisson::MultiPoisson()
  : PDFunction(),
    _N(vector<double>()),
    _Ngen(vector<double>()),
    _efl(vector<vector<double> >()),
    _bkg(vector<vector<double> >()),
    _meanefl(vector<double>()),
    _meanbkg(vector<double>()),    
    _random(TRandom3()),
    _nbins(0),
    _index(-1),
    _profile(false)
{}


MultiPoisson::MultiPoisson(string filename)
  : PDFunction(),
    _N(vector<double>()),
    _Ngen(vector<double>()),
    _efl(vector<vector<double> >()),
    _bkg(vector<vector<double> >()),
    _meanefl(vector<double>()),
    _meanbkg(vector<double>()),    
    _random(TRandom3()),
    _nbins(0),
    _index(-1),
    _profile(false) 
{
  // open input text file with format
  // bin1   bin2   ... 
  // count1 count2 ... 
  // efl1   efl2   ...
  // bkg1   bkg2   ...
  
  ifstream inp(filename);
  if ( ! inp.good() )
    {
      cout << "** MultiPoisson - unable to open file " << filename << endl;
      exit(0);
    }
  
  // read lines
  vector<string> records;
  string line;
  while(getline(inp, line))
    {
      if (line.substr(0,1) == "#") continue;
      records.push_back(line);
    }

  // get header and count number of columns
  vector<string> header;
  stringstream hin(records[0]);
  while (hin >> line)
    {
      header.push_back(line);
      //cout << line << endl;
    }
  _nbins = static_cast<int>(header.size());
  //cout << "number of bins = " << _nbins << endl;

  // get counts (record 1)
  vector<double> efl;
  vector<double> bkg;
  stringstream din(records[1]);
  double x;
  for(int i=0; i < _nbins; i++)
    {
      din >> x;
      _N.push_back(x);
      _meanefl.push_back(0);
      _meanbkg.push_back(0);
      efl.push_back(0);
      bkg.push_back(0);
      //cout << "count[" << i << "] = " << x << endl;
    }
  
  // loop over sample
  int sampleSize = (records.size()-2)/2;
  for(int c=0; c < sampleSize; c++)
    {
      int offset = 2*c + 1;
      
      // get effective luminosities
      stringstream ein(records[offset+1]);
      for(int i=0; i < _nbins; i++) ein >> efl[i];
     
      // get backgrounds
      stringstream bin(records[offset+2]);
      for(int i=0; i < _nbins; i++) bin >> bkg[i];

      add(efl, bkg);
    }
  computeMeans();
}

MultiPoisson::MultiPoisson(vector<double>& N)
  : PDFunction(),
    _N(N),
    _Ngen(N),
    _efl(vector<vector<double> >()),
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

 void MultiPoisson::add(vector<double>& efl, vector<double>& bkg)
{
  _efl.push_back(efl);
  _bkg.push_back(bkg);
}

void MultiPoisson::update(int ii, vector<double>& efl)
{
  if ( ii < 0 ) return;
  if ( ii > (int)(_efl.size()-1) ) return;
  copy(efl.begin(), efl.end(), _efl[ii].begin());
}

void MultiPoisson::computeMeans()
{
  int M = _efl.size();
  for(int ibin=0; ibin < _nbins; ++ibin)
    {
      _meanefl[ibin] = 0.0;
      _meanbkg[ibin] = 0.0;
      for(int ii=0; ii < M; ++ii)
	{
	  _meanefl[ibin] += _efl[ii][ibin];
	  _meanbkg[ibin] += _bkg[ii][ibin];
	}
      _meanefl[ibin] /= M;
      _meanbkg[ibin] /= M;
    }  
}

void MultiPoisson::set(int ii)
{
  if ( ii < 0 ) return;
  if ( ii > (int)(_efl.size()-1) ) return;
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
  if ( ii > (int)(_efl.size()-1) ) return vector<double>();
  vector<double> d(2*_efl[ii].size());
  double* a = &d[0];
  copy(_efl[ii].begin(), _efl[ii].end(), a);
  copy(_bkg[ii].begin(), _bkg[ii].end(), a + _efl[ii].size());
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
  int nconstants = _efl.size();
  int icon = _random.Integer(nconstants-1);      
  for(int ibin=0; ibin < _nbins; ++ibin)
    {
      double mean = _efl[icon][ibin] * sigma + _bkg[icon][ibin];
      _Ngen[ibin] = _random.Poisson(mean);
    }
  return _Ngen;
}

double 
MultiPoisson::operator() (std::vector<double>& N, double sigma)
{
  int first = 0;
  int last  = _efl.size()-1;
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
	      double mean = _efl[icon][ibin] * sigma + _bkg[icon][ibin];
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
