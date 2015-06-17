#ifndef PDFUNCTION_H
#define PDFUNCTION_H
//-------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: PDFunction.h
// Description: Base class for probability density functions, 
// modeled as function objects (that is, functors). This base 
// class is abstract and therefore cannot be instantiated. 
// Its purpose is
// 1) to enforce the same interface on all derived classes and
// 2) to allow for a generic algorithm, such as implemented in
//    the ReferencePrior class, to be applicable to any derived 
//    class, such as PoissonGamma. 
//
// Created: June 11, 2010
// Modifications: 
//
//--------------------------------------------------------------
#include <vector>
#include <iostream>
#include <stdlib.h>
#include "TROOT.h"

/**  Base class for probability density functions, modeled as 
     function objects (that is, functors). This base class 
     is abstract and therefore cannot be instantiated. Its 
     purpose is to enforce the same interface on all 
     derived classes.
*/
class PDFunction
{
 public:
  
  ///
  PDFunction() {}
  
  /** This forces the destructor of the derived class 
      to be called as well as this destructor
  */
  virtual ~PDFunction() {}
  
  /** The following pure virtual functions must be 
      overridden by derived classes. 
      See for example, MarginalizedPoissonGamma
  */
  virtual std::vector<double>& generate(double theta)=0;
  
  /** The following pure virtual functions must be 
      overridden by derived classes. 
      See for example, MarginalizedPoissonGamma
  */
  virtual double operator() (std::vector<double>& data, double theta)=0; 

 private:
  ClassDef(PDFunction,1)
};


#endif
