#ifndef ExpectedLimits_H
#define ExpectedLimits_H
//--------------------------------------------------------------
//
// File: ExpectedLimits.cc
// Description: Standalone expected limits calculator.
// 
// Created: 11 Jan 2011 Harrison B. Prosper
// Updated: 26 May 2017 HBP implement
//--------------------------------------------------------------
#include <vector>
#include <string>
#include "PDFunction.h"
#include "LimitCalculator.h"
//--------------------------------------------------------------
/** Compute expected limits.
    <p>
    Generate 
 */
class ExpectedLimits
{
public:
  static std::vector<double> dummy;
  
  ExpectedLimits ();

  /** Compute quantiles of limits distribution, generated internally.
      @param ensemble size - size of ensemble of limits
      @param calculator  - limit calculator (Bayes or Wald)
  */
  ExpectedLimits(LimitCalculator& calculator,
		 int ensemblesize=200,
		 std::vector<double>& prob_=dummy);

  virtual ~ExpectedLimits();

  virtual std::vector<double> prob() { return _prob; }
  
  virtual std::vector<double> operator() (double mu=1, int ensemblesize=-1);
  
private:
  LimitCalculator* _calculator;
  int _ensemblesize;
  std::vector<double> _prob;
  std::vector<double> _limit;

  int _debuglevel;
};

#endif


