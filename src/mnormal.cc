////////////////////////////////////////////////////////////////////////////
// File: mnormal.cc
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
// Updated: 17-Nov-2012 HBP add methods to return Gaussian variates and
//                          to re-use it
////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib> 
#include "TRandom3.h"
#include "mnormal.h"

using namespace std;

mnormal::mnormal() 
  : n(0) {}

TMatrixDSym
mnormal::covariance(TMinuit& minuit, int N)
{
  int ndim = N * N;
  double* errmat = new double(ndim);
  minuit.mnemat(errmat, ndim);

  TMatrixDSym cov(N);
  for(int ii=0; ii < N; ++ii)
    for(int jj=0; jj < N; ++jj)
      cov(ii, jj) = errmat[ii*ndim+jj];
  delete errmat;
  return cov;
}

mnormal::mnormal(std::vector<double>& ai)
    : random(TRandom3()),
      n(ai.size()),
      a(ai)
{
  v.clear();
  z.clear();
  c.clear();
  for (int i = 0; i < n; i++) z.push_back(0);
}

mnormal::mnormal(std::vector<double>& ai,
		 TMatrixDSym& cov)
  : random(TRandom3()),
    n(ai.size()),
    a(ai)
{
  v.clear();
  z.clear();
  c.clear();
  
  // store covariance matrix
  vector<double> row(n);
  for (int i = 0; i < n; i++)
    { 
      z.push_back(0);
      for (int j = 0; j < n; j++)
	row[j] = cov[i][j];
      addRow(row);
    }
}

void mnormal::setSeed(int seed)
{
  random.SetSeed(seed);
}

void mnormal::addRow(vector<double>& row)
{
  v.push_back(row);
  
  if ( v.size() < (unsigned int)n ) return;
  
  for (int i = 0; i < n; i++) c.push_back(z);

  // Compute Cholesky square root of covariance matrix
  // v = c*c^T
  ////////////////////////////////////////////////////
  for (int j = 0; j < n; j++)
    {
      // Compute diagonal terms
      /////////////////////////
      double y = 0;
      for (int k = 0; k < j; k++) y += c[j][k]*c[j][k];
      double x = v[j][j]-y;
      if      ( x <  0.0 ) 
	{
	  cout << "** ERROR ** Matrix not positive definite\n";
	  cout << "   Need to increase diagonal element "
               << "v(" << j+1 << "," << j+1 << ") = "
	       << v[j][j] 
	       << " by > " << fabs(x) << endl;
	  exit(0);
	}
      else if ( x == 0.0 )
	{
	  cout << "** ERROR ** Matrix singular\n";
	  cout << "   Need to add an offset " 
                 << "to diagonal element v(" << j+1 << "," << j+1 << ") = "
	       << v[j][j] << endl;
	  exit(0);
	}
      
      c[j][j] = sqrt(x);
      
      // Compute off-diagonal terms
      /////////////////////////////
      for (int i = j; i < n; i++)
	{
	  double yy = 0;
	  for (int k = 0; k < j; k++) yy += c[i][k]*c[j][k];
	  c[i][j] = (v[i][j] - yy)/c[j][j];
	}
    }
}

mnormal::~mnormal() {}

vector<double>& mnormal::getZ() { return z; }

// Generate pseudo-random vectors
/////////////////////////////////
bool mnormal::generate(std::vector<double>& x)
{
  bool positive = true;
  for (int i = 0; i < n; i++) z[i] = random.Gaus();
  
  for (int i = 0; i < n; i++)
    {
      double y = 0;
      for (int j = 0; j < n; j++) y += c[i][j]*z[j];
      y = y + a[i];
      x[i] = y;
      if ( x[i] < 0.0 ) positive = false;
    }
  return positive;
}

// Generate pseudo-random vectors
/////////////////////////////////
bool mnormal::generate(std::vector<double>& x, std::vector<double>& Z)
{
  bool positive = true;
  
  if ( Z.size() != z.size() )
    {
      cout << "Z size mismatch!" << endl;
      exit(0);
    }
  
  for (int i = 0; i < n; i++) z[i] = Z[i];
  
  for (int i = 0; i < n; i++)
    {
      double y = 0;
      for (int j = 0; j < n; j++) y += c[i][j]*z[j];
      y = y + a[i];
      x[i] = y;
      if ( x[i] < 0.0 ) positive = false;
    }
  return positive;
}


void mnormal::printme()
{
  cout << "\nmnormal: Cholesky Square Root of Matrix" << endl;
  
  char record[80];
  for (unsigned int i = 0; i < c.size(); i++)
    {
      for (unsigned int j = 0; j < c.size(); j++)
	{
	  sprintf(record, " %10.3e", c[i][j]);
	  cout << record;
	}
      cout << endl;
    }
  
  vector<float> zero(n);
  for (int i = 0; i < n; i++) zero[i] = 0;
  
  vector<vector<float> > b(n);
  for (int i = 0; i < n; i++) b[i] = zero;
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
	b[i][j] += c[i][k]*c[j][k];
  
  cout << "\nmnormal: Original Matrix\n";
  
  for (unsigned int i = 0; i < v.size(); i++)
    {
      for (unsigned int j = 0; j < v[i].size(); j++)
	{
	  sprintf(record, " %10.3e", v[i][j]);
	  cout << record; 
	}
      cout << endl;
    }
  
  cout << "\nmnormal: Reconstructed Matrix\n";
  
  for (unsigned int i = 0; i < b.size(); i++)
    {
      for (unsigned int j = 0; j < b[i].size(); j++)
	{
	  sprintf(record, " %10.3e", b[i][j]);
	  cout << record; 
	} 
      cout << endl;
    }
}
