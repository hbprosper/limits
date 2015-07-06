#ifndef MULTIPOISSON_H
#define MULTIPOISSON_H
//--------------------------------------------------------------
// File: MultiPoisson.h
// Description: Implement the multi-Poisson model averaged
//              with respect to a prior specified as a swarm
//              of points. Model assumes mean = efl * sigma + bkg
//              for each bin, where efl = eff * luminosity.
// 
// Created: 11-Jun-2010 Harrison B. Prosper & Supriya Jain
// Updated: 11-Aug-2014 HBP - add option to profile multi-Poisson
//                      model (not yet implemented!).
//          05-Jul-2015 HBP - assume user supplies the effective
//                      luminosity efl = eff * L rather than
//                      eff and L separately. This makes for a
//                      cleaner implementation.
//--------------------------------------------------------------
#include <vector>
#include <algorithm>
#include "TRandom3.h"
#include "PDFunction.h"

/** Implement the multi-Poisson model averaged over evidence-based prior.
*/
class MultiPoisson : public PDFunction
{
 public:
  ///
  MultiPoisson();

  /** Default constructor.  
      @param filename - name of text file containing counts, effective
      luminosities, and backgrounds.
  */
  MultiPoisson(std::string filename);
 
  
  /** Main constructor. 
      @param N - observed counts
  */
  MultiPoisson(std::vector<double>& N);
 
  ///
  ~MultiPoisson();
	     
  /** Generate data for one experiment.
      @param sigma - value of parameter of interest
  */
  std::vector<double>& generate(double sigma);
  
  /** Compute likelihood.
      @param N - observed data
      @param sigma - value of parameter of interest 
  */
  double operator() (std::vector<double>& N, double sigma);

  /** Compute likelihood using internally cached data.
      @param sigma - value of parameter of interest 
  */
  double operator() (double sigma);

  /** If true, profile rather than average.
      Not yet implemented.
   */
  void profile(bool yes=true) {_profile=yes;}

  /** Add one set of efficiency, background, and luminosity parameters.
      @param efl - effective luminosity parameters (efl = eff * lumi)
      @param bkg - background parameters
   */
  void add(std::vector<double>& eff, std::vector<double>& bkg);

  /** Update effective luminosities.
   */
  void update(int ii, std::vector<double>& efl);
  
  /** Compute mean effective luminosities and backgrounds.
   */
  void computeMeans();
  
  void set(int ii);
  void reset();
  void setSeed(int seed);

  ///
  std::vector<double> get(int ii);

  ///
  std::vector<double> counts() { return _N; }

  /// Return average effective luminosities (efficiency * luminosity).
  std::vector<double> eluminosity() { return _meanefl; }

  /// Return average backgrounds.
  std::vector<double> background() { return _meanbkg; }

  /// Return sample size.
  int size() { return _efl.size(); }
  
 private:
    std::vector<double> _N;
    std::vector<double> _Ngen;
    std::vector<std::vector<double> > _efl;
    std::vector<std::vector<double> > _bkg;
    std::vector<double> _meanefl;
    std::vector<double> _meanbkg;
    
    TRandom3 _random;
    int _nbins;
    int _index;
    bool _profile;
};

#endif

