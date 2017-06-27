//--------------------------------------------------------------
// File: MultiPoissonGamma.cc
// Description: Implement the multi-Poisson-Gamma model
//              averaged over an evidence-based prior.
// 
// Created: 11-Jun-2014 HBP
// Updated  12-Jun-2016 HBP allow input of histograms
//          25-May-2017 HBP simplify config file format
//                      warning: not backwards compatible!
//--------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TError.h"
#include "MultiPoissonGamma.h"

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
  // examine file and determine whether we have histogram input)
  ifstream inp(filename.c_str());
  if ( ! inp.good() )
    {
      cout << "** MultiPoissonGamma - unable to open file "
	   << filename << endl;
      exit(0);
    }
  
  // read lines
  vector<string> records;
  string line;
  bool useHistograms = false;
  
  while(getline(inp, line))
    {
      // remove blank lines
      line = strip(line);
      if (line == "") continue;
      
      // remove lines that start with #
      if (line.substr(0,1) == "#") continue;

      // check if a line contains the name of a root file
      // if so, assume histogram input has been specified for
      // all input data
      TString str(line.c_str());
      if ( str.Contains(".root") ) useHistograms = true;
      
      records.push_back(line);
    }

  if (records.size() == 0)
    {
      Error("MultiPoisson",
	    "zero uncommented records in file %s", filename.c_str());
      exit(0);
    }
    
  // now decode 
  if ( useHistograms )
    _readRootFile(records);
  else
    _readTextFile(records);
}

void
MultiPoissonGamma::_readTextFile(vector<string>& records)
{
  // open input text file with format
  // open input text file with format
  // number of bins  ... 
  // count1 count2 ...
  // number of points
  // S1   S2   ...
  // dS1  dS2  ...
  // B1   B2   ...
  // dB1  dB2  ...

  // get number of bins
  _nbins = 0;
  int lineno = 0;
  stringstream hin(records[lineno]);
  try
    {
      hin >> _nbins;
    }
  catch (...)
    {
      Error("MultiPoisson",
	    "unable to get number of bins; check file format");
      exit(0);
    }  

  // get observed counts
  vector<double> sig;
  vector<double> dsig;
  vector<double> bkg;
  vector<double> dbkg;
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
	  Error("MultiPoissonGamma",
		"problem accessing observed count at bin %d", c);
	  exit(0);	  
	}
      _N.push_back(x);
      _Ngen.push_back(0);
      sig.push_back(0);
      bkg.push_back(0);
      dsig.push_back(0);
      dbkg.push_back(0);
    }
  
  // get number of samples
  int samplesize = 0;
  lineno++;
  stringstream nin(records[lineno]);
  try
    {
      nin >> samplesize;
    }
  catch (...)
    {
      Error("MultiPoissonGamma",
	    "unable to get sample size; check file format");
      exit(0);
    }  

  // loop over sampled points  
  for(int ii=0; ii < samplesize; ii++)
    {
      try
	{
	  // get signals
	  lineno++;
	  stringstream sin(records[lineno]);
	  for(int i=0; i < _nbins; i++) sin >> sig[i];
      
	  // get associated uncertainties
	  lineno++;
	  stringstream dsin(records[lineno]);
	  for(int i=0; i < _nbins; i++) dsin >> dsig[i];
	}
      catch (...)
	{
	  Error("MultiPoisson",
		"problem decoding signal counts:\n%s",
		records[lineno].c_str());
	  exit(0);      
	}

      try
	{
	  // get backgrounds
	  lineno++;
	  stringstream bin(records[lineno]);
	  for(int i=0; i < _nbins; i++) bin >> bkg[i];
      
	  // get background uncertainties
	  lineno++;
	  stringstream dbin(records[lineno]);
	  for(int i=0; i < _nbins; i++) dbin >> dbkg[i];
	}
      catch (...)
	{
	  Error("MultiPoisson",
		"problem decoding background counts:\n%s",
		records[lineno].c_str());
	  exit(0);	  
	}
      add(sig, dsig, bkg, dbkg);
    }
}

void
MultiPoissonGamma::_readRootFile(vector<string>& records)
{
  // open input text file with format
  // number of bins
  // root data filename      histogram name
  // number of sampled points
  // root signal filename    histogram name  [rel-uncertainty]
  // root backgroud filename histogram name  [rel-uncertainty] 
  //    :  :
  // get number of bins
  _nbins = 0;
  int lineno = 0;
  stringstream hin(records[lineno]);
  try
    {
      hin >> _nbins;
    }
  catch (...)
    {
      Error("MultiPoissonGamma",
	    "unable to get number of bins; check file format");
      exit(0);
    }

  // Get observed counts
  // rfilename = root file name
  // histname  = histogram name
  lineno++;
  istringstream iin(records[lineno]);
  string rfilename, histname;
  iin >> rfilename >> histname;

  vector<double> c;
  vector<double> dc;
  _getcounts(rfilename, histname, c, dc);
  if (_nbins != (int)c.size())
    Warning("MultiPoissonGamma",
	    "number of bins specified %d != number of bins %d "
	    "in data histogram\n"
	    "will use smaller bin count", _nbins, (int)c.size());
  _nbins = min(_nbins, (int)c.size());
  
  // open output text file
  string ofilename("limits.dat");
  ofstream fout(ofilename.c_str());
  
  // write out number of bins
  fout << "# automatically generated by MultiPoissonGamma" << endl;
  fout << "# number of bins" << endl;
  fout << _nbins << endl;
  
  char arecord[10000];
  
  // write out observed counts
  fout << "# observed counts" << endl;
  for(int i=0; i < _nbins; i++)
    {
      sprintf(arecord, " %9.0f", c[i]);
      fout << arecord;
    }
  fout << endl;
  
  // get number of samples
  int samplesize = 0;
  lineno++;
  stringstream nin(records[lineno]);
  try
    {
      nin >> samplesize;
    }
  catch (...)
    {
      Error("MultiPoissonGamma",
	    "unable to get sample size; check file format");
      exit(0);
    }
  
  fout << "# number of sampled points" << endl;
  fout << samplesize << endl;
  
  // loop over sample
  
  for(int ii=0; ii < samplesize; ii++)
    {
      // write out signals      
      lineno++;
      istringstream sin(records[lineno]);
      sin >> rfilename >> histname;
      _getcounts(rfilename, histname, c, dc);
      if ( (int)c.size() < _nbins )
	{
	  Error("MultiPoissonGamma",
		"bin count %d for histogram %s/%s < %d",
		(int)c.size(), 
		rfilename.c_str(),
		histname.c_str(),
		_nbins);
	  exit(0);
	}
      
      fout << "# " << ii+1 <<  " signal " << endl;
      
      // write out counts
	
      for(int i=0; i < _nbins; i++)
	{
	  sprintf(arecord, " %9.3e", c[i]);
	  fout << arecord;
	}
      fout << endl;
      fout << "# signal uncertainties" << endl;
      
      // write out uncertainties
      for(int i=0; i < _nbins; i++)
	{
	  sprintf(arecord, " %9.3e", dc[i]);
	  fout << arecord;
	}
      fout << endl;

      // write out backgrounds
      lineno++;
      istringstream bin(records[lineno]);
      bin >> rfilename >> histname;
      _getcounts(rfilename, histname, c, dc);
      if ( (int)c.size() != _nbins )
	{
	  Error("MultiPoissonGamma",
		"bin count %d for histogram %s/%s < %d",
		(int)c.size(), 
		rfilename.c_str(),
		histname.c_str(),
		_nbins);
	  exit(0);
	}

      fout << "# " << ii+1 << " backgrounds" << endl;
      
      // write out counts
      for(int i=0; i < _nbins; i++)
	{
	  sprintf(arecord, " %9.3e", c[i]);
	  fout << arecord;
	}
      fout << endl;
      fout << "# background uncertainties" << endl;
      
      // write out uncertainties
      for(int i=0; i < _nbins; i++)
	{
	  sprintf(arecord, " %9.3e", dc[i]);
	  fout << arecord;
	}
      fout << endl;      
    }

  // -------------------------------------------------------  
  // now read text file that's just been created
  // -------------------------------------------------------
  ifstream inp(ofilename.c_str());
  if ( ! inp.good() )
    {
      Error("MultiPoissonGamma",
	    "unable to open file %s", ofilename.c_str());
      exit(0);
    }
  
  vector<string> orecords;
  string line;
  while(getline(inp, line))
    {
      // remove blank lines
      line = strip(line);
      if (line == "") continue;
      // remove lines that start with #
      if (line.substr(0,1) == "#") continue;
      orecords.push_back(line);
    }
  _readTextFile(orecords);  
}

void
MultiPoissonGamma::_getcounts(string rfilename, string histname,
			      vector<double>& c, vector<double>& dc)
{
  c.clear();
  dc.clear();
  
  TFile rfile(rfilename.c_str());
  if ( ! rfile.IsOpen() )
    {
      Error("MultiPoissonGamma", "unable to open file %s",
	    rfilename.c_str());
      exit(0);
    }
  TH1* h = (TH1*)rfile.Get(histname.c_str());
  if ( ! h )
    {
      Error("MultiPoissonGamma",
	    "unable to get histogram %s from file %s",
	    histname.c_str(), rfilename.c_str());
      exit(0);
    }

  if ( h->InheritsFrom("TH2") )
    {
      for(int i=0; i < h->GetNbinsX(); i++)
	for(int j=0; j < h->GetNbinsY(); j++)
	  {
	    c.push_back(h->GetBinContent(i+1, j+1));
	    dc.push_back(h->GetBinError(i+1, j+1));
	  }
    }
  else
    {
      for(int i=0; i < h->GetNbinsX(); i++)
	{
	  c.push_back(h->GetBinContent(i+1));
	  dc.push_back(h->GetBinError(i+1));
	}
    }
  rfile.Close();
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

void MultiPoissonGamma::add(vector<double>& sig, vector<double>& dsig,
			    vector<double>& bkg, vector<double>& dbkg)
			    
			    
{
  vector<double> x;
  vector<double> a;
  _convert(sig, dsig, x, a);
  
  vector<double> y;
  vector<double> b;
  _convert(bkg, dbkg, y, b);

  _model.push_back( MultiPoissonGammaModel(_N, x, a, y, b) );
  if ( _model.size() % 50 == 0 )
    cout << "=> MultiPoissonGamma: added "
	 << _model.size() << " distributions" << endl;
}

void MultiPoissonGamma::update(int ii,
			       vector<double>& sig,
			       vector<double>& dsig)
{
  if ( ii < 0 ) return;
  if ( ii > (int)(_model.size()-1) ) return;

  vector<double> x;
  vector<double> a;
  _convert(sig, dsig, x, a);
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
MultiPoissonGamma::generate(double mu)
{
  if(_nbins == 0)
    {
      cout << "MultiPoissonGamma::generate: nbins = 0" << endl;
      exit(0);
    }

  int ii = _random.Integer(_model.size()-1);
  _Ngen  = _model[ii].generate(mu);
  return _Ngen;
}


double 
MultiPoissonGamma::operator() (std::vector<double>& N, double mu)
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
	  double p = _model[ii](N, mu);
	  likelihood += p;
	}
      likelihood /= nconstants;
    }
  return likelihood;
}

void 
MultiPoissonGamma::setSeed(int seed) { _random.SetSeed(seed); }

double 
MultiPoissonGamma::operator() (double mu)
{
  return (*this)(_N, mu);
}

void MultiPoissonGamma::_convert(vector<double>& sig, vector<double>& dsig,
				 vector<double>& x, vector<double>& a)
{
  x.clear();
  a.clear();
  for(size_t i=0; i < sig.size(); i++)
    {
      double c  = sig[i];
      double dc = dsig[i];
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
