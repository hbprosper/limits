//--------------------------------------------------------------
// File: MultiPoissonGamma.cc
// Description: Implement the multi-Poisson-Gamma model
//              averaged over an evidence-based prior.
// 
// Created: June 11, 2014 HBP
// UpdatedL June 12, 2016 HBP allow input of histograms
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
      TString str(line.c_str());
      if ( str.Contains(".root") ) useHistograms = true;
      
      records.push_back(line);
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
  // bin1   bin2   ...
  // count1 count2 ...
  // sig1   sig2   ...
  // dsig1  dsig2   ...
  // bkg1   bkg2   ...
  // dbkg1  dbkg2   ...
 
  // get header and count number of columns
  vector<string> header;
  stringstream hin(records[0]);
  string line;
  while (hin >> line) header.push_back(line);
  int nbins = (int)header.size();
  
  // get counts (record 1)
  vector<double> sig;
  vector<double> dsig;
  vector<double> bkg;
  vector<double> dbkg;
  stringstream din(records[1]);
  double x;
  for(int i=0; i < nbins; i++)
    {
      din >> x;
      _N.push_back(x);
      sig.push_back(0);
      bkg.push_back(0);
      dsig.push_back(0);
      dbkg.push_back(0);
    }
  
  // loop over sample
  int sampleSize = (records.size()-2)/4;
  
  for(int c=0; c < sampleSize; c++)
    {
      int offset = 1 + 4*c;
      // get signals
      stringstream ein(records[offset+1]);
      for(int i=0; i < nbins; i++) ein >> sig[i];

      // get associated uncertainties
      stringstream dein(records[offset+2]);
      for(int i=0; i < nbins; i++) dein >> dsig[i];
 
      // get backgrounds
      stringstream bin(records[offset+3]);
      for(int i=0; i < nbins; i++) bin >> bkg[i];
      
      // get background uncertainties
      stringstream dbin(records[offset+4]);
      for(int i=0; i < nbins; i++) dbin >> dbkg[i];
 
      add(sig, dsig, bkg, dbkg);
    }
}

void
MultiPoissonGamma::_readRootFile(vector<string>& records)
{
  // open input text file with format
  // root data filename      histogram name
  // root signal filename    histogram name
  // root backgroud filename histogram name
  //    :  :

  // Get observed counts
  istringstream iin(records[0]);

  // rfilename = root file name
  // histname  = histogram name
  string rfilename, histname;
  iin >> rfilename >> histname;

  // get observed counst
  vector<double> c;
  vector<double> dc;
  _getcounts(rfilename, histname, c, dc);
  int nbins = (int)c.size();

  // open output text file
  string ofilename("limits.dat");
  ofstream fout(ofilename.c_str());
  
  // write out header
  char record[80];
  for(int i=0; i < nbins; i++)
    {
      sprintf(record, " bin%6.6d", i+1);
      fout << record;
    }
  fout << endl;
  
  // write out observed counts
  fout << "# observed counts" << endl;
  for(int i=0; i < nbins; i++)
    {
      sprintf(record, " %9.2f", c[i]);
      fout << record;
    }
  fout << endl;
    
  // loop over sample
  int sampleSize = (records.size()-1)/2;
  for(int ii=0; ii < sampleSize; ii++)
    {
      // write out signals
      istringstream in(records[ii+1]);
      in >> rfilename >> histname;
      _getcounts(rfilename, histname, c, dc);
      if ( (int)c.size() != nbins )
	{
	  cout << "** MultiPoissonGamma - mismatch in bin count" << endl
	       << "** for histogram " << rfilename << "/" << histname << endl;
	  exit(0);
	}
      
      fout << "# signals ====================== " << ii+1 << endl;
      
      // write out counts
      for(int i=0; i < nbins; i++)
	{
	  sprintf(record, " %9.2f", c[i]);
	  fout << record;
	}
      fout << endl;
      // write out uncertainties
      for(int i=0; i < nbins; i++)
	{
	  sprintf(record, " %9.2f", dc[i]);
	  fout << record;
	}
      fout << endl;

      fout << "# backgrounds" << endl;
      
      // write out backgrounds
      istringstream din(records[ii+2]);
      din >> rfilename >> histname;
      _getcounts(rfilename, histname, c, dc);
      if ( (int)c.size() != nbins )
	{
	  cout << "** MultiPoissonGamma - mismatch in bin count" << endl
	       << "** for histogram " << rfilename << "/" << histname << endl;
	  exit(0);
	}
      // write out counts
      for(int i=0; i < nbins; i++)
	{
	  sprintf(record, " %9.2f", c[i]);
	  fout << record;
	}
      fout << endl;
      // write out uncertainties
      for(int i=0; i < nbins; i++)
	{
	  sprintf(record, " %9.2f", dc[i]);
	  fout << record;
	}
      fout << endl;      
    }

  // -------------------------------------------------------  
  // now read text file that's just been created
  // -------------------------------------------------------
  ifstream inp(ofilename.c_str());
  if ( ! inp.good() )
    {
      cout << "** MultiPoissonGamma - unable to open file "
	   << ofilename << endl;
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
      cout << "** MultiPoissonGamma - unable to open file "
	   << rfilename << endl;
      exit(0);
    }
  TH1* h = (TH1*)rfile.Get(histname.c_str());
  if ( ! h )
    {
      cout << "** MultiPoissonGamma - unable to get histogram "
	   << histname << endl
	   << "** from file "
	   << rfilename << endl;
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
