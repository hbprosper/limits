//-------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
//         Luc Demortier       (luc@fnal.gov)
//         Supriya Jain        (sjain@fnal.gov)
// 
// 
// File: PDFunction.cc
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
#include "PDFunction.h"
ClassImp(PDFunction)

