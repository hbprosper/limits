#ifndef MNORMAL_H
#define MNORMAL_H
////////////////////////////////////////////////////////////////////////////
// File: mnormal.h
// Description: Generate a vector of variates according to a multi-variate
//              Gaussian.
// Usage:
//       (a) Initialization
//
//           mnormal r(a)
//                   Inputs:
//                      vector<double>         a   vector of mean values
//           for(unsigned int i=0; i < a.size(); i++)
//             {
//                      :   :
//               row[0] = ...
//
//               row[a.size()-1] = ...
//               r.addRow(row);   // Add ith row of covariance matrix 
//                      :   :
//             }
//
//       (b) Generation
//
//           ok = r.generate(x)
//                    Outputs:
//                      vector<double>         x   random vector
//                      bool                   ok  true if point is
//                                                 in positive "quadrant"
// Created: 6-Jun-2000 Harrison B. Prosper
//                     C++ version of my 1986 version of the routine!  
//
// Updated: 17-Nov-2012 HBP add methods to return Gaussian variates so that
//                          we can re-use them
////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib> 
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMinuit.h"

namespace {
  typedef std::vector<std::vector<double> > vvdouble; 
  std::vector<double> ZG;
};

/// Generate variates according to a multivariate Gaussian.
class mnormal
{
public:
  ///
  mnormal();
  
  /// Constructor with vector of means.
  mnormal(std::vector<double>& ai);

  /// Vector of means and covariance matrix.
  mnormal(std::vector<double>& ai,
	  TMatrixD& cov);
  
  /// Vector of means, vector of standard deviations and correlation matrix.
  mnormal(std::vector<double>& ai,
	  std::vector<double>& stddev,
	  TMatrixD& cor);
  ///
  void setSeed(int seed);

  /// Add one row of covariance matrix.
  void addRow(std::vector<double>& row);

  ~mnormal();

  /// Get NxN covariance matrix from Minuit.
  static TMatrixD covariance(TMinuit& minuit, int N);
  
  /// Get vector of unit variance, zero mean, variates.
  std::vector<double>& getZ();

  /// Generate random vectors of Gaussian variates.
  bool generate(std::vector<double>& x, std::vector<double>& Z=ZG);

  /// Print Cholesky square root of covariance matrix.
  void printme();

private:
  TRandom3 random;
  int n;
  std::vector<double> a;
  std::vector<double> z;
  std::vector<std::vector<double> > v;
  std::vector<std::vector<double> > c;
  
};
#endif
