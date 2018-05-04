//--------------------------------------------------------------
// File: ExpectedLimits.cc
// Description: Standalone expected limits calculator.
// 
// Created: 11 Jan 2011 Harrison B. Prosper
// Updated: 26 May 2017 HBP implement
//--------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include "ExpectedLimits.h"

using namespace std;
// ---------------------------------------------------------------------------
vector<double> ExpectedLimits::dummy;

ExpectedLimits::ExpectedLimits()
  : _calculator(0),
    _ensemblesize(0),
    _prob(dummy),
    _debuglevel(0)
{}


ExpectedLimits::ExpectedLimits(LimitCalculator& calculator,
			       int ensemblesize,
			       std::vector<double>& prob)
			       
  : _calculator(&calculator),
    _ensemblesize(ensemblesize),
    _prob(prob),
    _limit(vector<double>(ensemblesize)),
    _rms(0),
    _debuglevel(0)
{
  if ( getenv("DBExpectedLimits") > 0 )
    _debuglevel = atoi(getenv("DBExpectedLimits"));
  else
    _debuglevel = 0;
  
  if (_prob.size() == 0)
    {
      _prob.push_back(0.0230);
      _prob.push_back(0.1579);
      _prob.push_back(0.5000);
      _prob.push_back(0.8415);
      _prob.push_back(0.9770);
    }
}


ExpectedLimits::~ExpectedLimits()
{
}

vector<double>
ExpectedLimits::operator()(double true_value, bool compute_rms)
{
  char record[80];
  int step = _ensemblesize / 4;
  _rms = 0;
  for(int c=0; c < _ensemblesize; c++)
    {
      if ( c % step == 0 ) cout << "\tgenerating sample:\t" << c;
      
      // generate a data set assuming the background only hypothesis
      // that is, mu=0
      vector<double>& d = _calculator->pdf()->generate(true_value);
      if ( _debuglevel > 2 )
	{
	  cout << endl << c << "\tgenerated data: " << endl;
	  for(size_t ii=0; ii < d.size(); ii++)
	    {
	      sprintf(record, " %9.0f", d[ii]);
	      cout << record;
	    }
	  cout << endl;
	}
      
      // update data in calculator
      _calculator->setData(d);
      
      // compute 95% limit
      _limit[c] = _calculator->percentile();

      if ( compute_rms )
	{
	  double estimate = _calculator->estimate();
	  double de = (true_value - estimate);
	  _rms += de*de;
	}
      if ( c % step == 0 )
	{
	  if ( compute_rms )
	    {
	      if ( c > 0 )
		cout << "\trms = " << sqrt(_rms / c) << endl;
	      else
		cout << endl;
	    }
	  else
	    cout << endl;
	}
    }

  if ( compute_rms )
    _rms = sqrt(_rms / _ensemblesize);
		
  // now sort limits in increasing order
  sort(_limit.begin(), _limit.end());

  // get percentiles
  vector<double> percentiles(_prob.size());
  for(size_t ii=0; ii < _prob.size(); ii++)
    {
      // compute ordinal value of percentile
      double q = _prob[ii] * _ensemblesize;
      int    i = (int)q;
      double x = q - i;
      percentiles[ii] = x * _limit[i+1] + (1 - x) * _limit[i]; 
    }
  return percentiles;
}
