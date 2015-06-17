#ifndef PRIORFUNCTION_H
#define PRIORFUNCTION_H
//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: PriorFunction.h
// Description: Base class for prior density functions, 
//              modeled as function objects.
//              This base class is abstract and therefore 
//              cannot be instantiated. Its purpose is to enforce
//              the same interface on all derived classes.
//
// Created: June 11, 2010
// Modifications: 
//
//--------------------------------------------------------------
#include <vector>
#include <iostream>
#include <stdlib.h>

/**  Base class for prior density functions, modeled as 
     function objects. This base class 
     is abstract and therefore cannot be instantiated. Its 
     purpose is to enforce the same interface on all 
     derived classes.
*/
class PriorFunction
{
 public:
  
  ///
  PriorFunction() {}
	
  /** This forces the destructor of the derived class 
      to be called as well as this destructor
  */
  virtual ~PriorFunction() {}

  /** The following pure virtual functions must be 
      overridden by derived classes. 
      See for example, ReferencePrior
  */
  virtual double operator() (double poi)=0; 
};


#endif

