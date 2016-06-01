//--------------------------------------------------------------
// File: MultiPoissonGamma.cc
// Description: Implement the multi-Poisson-Gamma model
//              averaged over an evidence-based prior.
// 
// Created: June 11, 2014
//--------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "TMath.h"
#include "MultiPoissonGamma.h"

using namespace std;

MultiPoissonGamma::MultiPoissonGamma()
  : PDFunction(),
    _N(vector<double>()),
    _Ngen(vector<double>()),
    _model(vector<MultiPoissonGammaModel>()),
    _random(TRandom3()),
    _nbins(0),
    _index(-1),
    _profile(false)
{}


MultiPoissonGamma::MultiPoissonGamma(string filename)
  : PDFunction(),
    _N(vector<double>()),
    _Ngen(vector<double>()),
    _model(vector<MultiPoissonGammaModel>()),
    _random(TRandom3()),
    _nbins(0),
    _index(-1),
    _profile(false) 
{
  // open input text file with format
  // bin1   bin2   ...
  // count1 count2 ...
  // efl1   efl2   ...
  // defl1  defl2   ...
  // bkg1   bkg2   ...
  // dbkg1  dbkg2   ...
  
  ifstream inp(filename.c_str());
  if ( ! inp.good() )
    {
      cout << "** MultiPoissonGamma - unable to open file " << filename << endl;
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
  while (hin >> line) header.push_back(line);
  int nbins = (int)header.size();
  
  // get counts (record 1)
  vector<double> efl;
  vector<double> defl;
  vector<double> bkg;
  vector<double> dbkg;
  stringstream din(records[1]);
  double x;
  for(int i=0; i < nbins; i++)
    {
      din >> x;
      _N.push_back(x);
      efl.push_back(0);
      bkg.push_back(0);
      defl.push_back(0);
      dbkg.push_back(0);
    }
  
  // loop over sample
  int sampleSize = (records.size()-2)/4;
  
  for(int c=0; c < sampleSize; c++)
    {
      int offset = 1 + 4*c;
      // get effective luminosities
      stringstream ein(records[offset+1]);
      for(int i=0; i < nbins; i++) ein >> efl[i];

      // get effective luminosity uncertainties
      stringstream dein(records[offset+2]);
      for(int i=0; i < nbins; i++) dein >> defl[i];
 
      // get backgrounds
      stringstream bin(records[offset+3]);
      for(int i=0; i < nbins; i++) bin >> bkg[i];
      
      // get background uncertainties
      stringstream dbin(records[offset+4]);
      for(int i=0; i < nbins; i++) dbin >> dbkg[i];
 
      add(efl, defl, bkg, dbkg);
    }
}

MultiPoissonGamma::MultiPoissonGamma(vector<double>& N)
  : PDFunction(),
    _N(N),
    _Ngen(N),
    _random(TRandom3()),
    _nbins((int)N.size()),
    _index(-1),
    _profile(false)
{}

MultiPoissonGamma::~MultiPoissonGamma() 
{}

void MultiPoissonGamma::add(vector<double>& efl, vector<double>& defl,
			    vector<double>& bkg, vector<double>& dbkg)
			    
			    
{
  vector<double> x;
  vector<double> a;
  _convert(efl, defl, x, a);
  
  vector<double> y;
  vector<double> b;
  _convert(bkg, dbkg, y, b);

  _model.push_back( MultiPoissonGammaModel(_N, x, a, y, b) );
}

void MultiPoissonGamma::update(int ii,
			       vector<double>& efl,
			       vector<double>& defl)
{
  if ( ii < 0 ) return;
  if ( ii > (int)(_model.size()-1) ) return;

  vector<double> x;
  vector<double> a;
  _convert(efl, defl, x, a);
  _model[ii].setX(x, a);
}

void MultiPoissonGamma::set(int ii)
{
  if ( ii < 0 ) return;
  if ( ii > (int)(_model.size()-1) ) return;
  _index = ii;
}

void MultiPoissonGamma::reset()
{
  _index = -1;
}

vector<double>&  
MultiPoissonGamma::generate(double sigma)
{
  if(_nbins == 0)
    {
      cout << "MultiPoissonGamma::generate: nbins = 0" << endl;
      exit(0);
    }

  int ii = _random.Integer(_model.size()-1);
  _Ngen  = _model[ii].generate(sigma);
  return _Ngen;
}


double 
MultiPoissonGamma::operator() (std::vector<double>& N, double sigma)
{
  int first = 0;
  int last  = _model.size()-1;
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
      for(int ii=first; ii <= last; ++ii)
	{
	  double p = _model[ii](N, sigma);
	  likelihood += p;
	}
      likelihood /= nconstants;
    }
  return likelihood;
}

void 
MultiPoissonGamma::setSeed(int seed) { _random.SetSeed(seed); }

double 
MultiPoissonGamma::operator() (double sigma)
{
  return (*this)(_N, sigma);
}

void MultiPoissonGamma::_convert(vector<double>& efl, vector<double>& defl,
				 vector<double>& x, vector<double>& a)
{
  x.clear();
  a.clear();
  for(size_t i=0; i < efl.size(); i++)
    {
      double c  = efl[i];
      double dc = defl[i];
      if ( c <= 0 ) c = 1.e-3;
      if ( dc<= 0 )dc = 1.e-4;
      double k = c / dc;
      k *= k;
      double k2 = k+2;
      double gamma = (k2 + sqrt(k2*k2 - 4))/2 ;
      double beta  = (sqrt(c*c + 4*dc*dc) - c)/2;
      x.push_back(gamma - 1);
      a.push_back(1.0/beta);
    }
}
