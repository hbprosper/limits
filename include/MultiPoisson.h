#ifndef MULTIPOISSON_H
#define MULTIPOISSON_H
//--------------------------------------------------------------
// File: MultiPoisson.h
// Description: Implement the multi-Poisson model averaged
//              with respect to a prior specified as a swarm
//              of points. Model assumes mean = eff*L + bkg for
//              each bin.
// 
// Created: 11-Jun-2010 Harrison B. Prosper & Supriya Jain
// Updated: 11-Aug-2014 HBP - add option to profile multi-Poisson
//                      model (not yet implemented!).
//--------------------------------------------------------------
#include <vector>
#include <algorithm>
#include "TRandom3.h"
#include "PDFunction.h"

/** Implements the multi-Poisson model
*/
class MultiPoisson : public PDFunction
{
 public:
  ///
  MultiPoisson();

  /** Constructor:  
      @param filename - name of text file containing counts, etc.
  */
  MultiPoisson(std::string filename, double L=1);
 
  
  /** Constructor:  
      @param N - observed counts
  */
  MultiPoisson(std::vector<double>& N);
 
  ///
  ~MultiPoisson();
	     
  /** Generate data for one experiment.
      @param sigma - value of parameter of interest
      @param useMean - if true, set effective luminosity and background to
      mean values.
  */
  std::vector<double>& generate(double sigma);

  /** Generate data for one experiment using the mean efficiency and background.
      @param sigma - value of parameter of interest
      @param useMean - if true, set effective luminosity and background to
      mean values.
  */
  std::vector<double>& generateUsingMeans(double sigma);

  
  /** Computes PDF.
      @param N - observed data
      @param sigma - value of parameter of interest 
  */
  double operator() (std::vector<double>& N, double sigma);

  double operator() (double sigma);

  /** if true, profile rather than average
   */
  void profile(bool yes=true) {_profile=yes;}

  void add(std::vector<double>& eff,
	   std::vector<double>& bkg,
	   double L=1);
  
  void update(int ii, std::vector<double>& eff);
  void computeMeans();
  void set(int ii);
  void reset();
  void setSeed(int seed);

  std::vector<double> get(int ii);

  ///
  std::vector<double> counts() { return _N; }

  /// Average efficiency X luminosity.
  std::vector<double> eluminosity() { return _meanefl; }

  /// Average background.
  std::vector<double> background() { return _meanbkg; }

  /// Sample size.
  int size() { return _L.size(); }
  
 private:
    std::vector<double> _N;
    std::vector<double> _Ngen;
    std::vector<double> _L;
    std::vector<std::vector<double> > _eff;
    std::vector<std::vector<double> > _bkg;

    std::vector<double> _meanefl;
    std::vector<double> _meanbkg;
    
    TRandom3 _random;
    int _nbins;
    int _index;
    bool _profile;
};

#endif

