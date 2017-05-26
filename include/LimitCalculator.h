#ifndef LimitCalculator_H
#define LimitCalculator_H
//--------------------------------------------------------------
//
// File: LimitCalculator.cc
// Description: Base class for limit calculators
// 
// Created: 11 Jan 2011 Harrison B. Prosper
// Updated: 26 May 2017 HBP define
//--------------------------------------------------------------
#include <vector>
#include <string>
#include "PDFunction.h"
#include "PriorFunction.h"
//--------------------------------------------------------------
/** Base class for limit calculators
 */
class LimitCalculator
{
public:
  ///
  LimitCalculator () {}

  ///
  virtual ~LimitCalculator () {}

  /// Return pdf of model.
  virtual PDFunction* pdf()=0;
  
  /// Compute percenttile of density.
  virtual double percentile(double p=-1)=0;

  ///
  virtual void setData(std::vector<double>& d)=0;
  
  /// Compute a Z-value
  virtual double zvalue(double mu=1)=0 ;
};

#endif


