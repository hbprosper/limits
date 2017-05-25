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
//          25-May-2017 HBP - use S and B instead of efl and bkg!
//--------------------------------------------------------------
#include <vector>
#include <algorithm>
#include "TRandom3.h"
#include "PDFunction.h"

/** Implement the multi-Poisson model averaged over an evidence-based prior.
    The evidence-based prior is given as a swarm of points over signal and
    backgrounds.
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
  virtual ~MultiPoisson();
	     
  /** Generate data for one experiment.
      @param mu - parameter of interest
  */
  std::vector<double>& generate(double mu);
  
  /** Compute likelihood.
      @param N - observed data
      @param mu - parameter of interest 
  */
  double operator() (std::vector<double>& N, double mu);

  /** Compute likelihood using internally cached data.
      @param  mu - value of parameter of interest 
  */
  double operator() (double mu);

  /** If true, profile rather than average.
      Not yet implemented.
   */
  void profile(bool yes=true) {_profile=yes;}

  /** Add one set of signal and background parameters.
      @param S - signals or effective luminosities (eff * lumi)
      @param B - backgrounds
   */
  void add(std::vector<double>& S, std::vector<double>& B);

  /** Update specified signal parameter point.
   */
  void update(int ii, std::vector<double>& S);
  
  /** Compute mean signals and backgrounds.
   */
  void computeMeans();
  
  void set(int ii);
  void reset();
  void setSeed(int seed);

  ///
  //std::vector<std::pair<double, double> > get(int ii);

  ///
  std::vector<double> counts() { return _N; }

  /// Return average signal or (efficiency * luminosity).
  std::vector<double> signal() { return _meanS; }

  /// Return average background.
  std::vector<double> background() { return _meanB; }

  /// Return sample size.
  int size() { return _S.size(); }
  
 private:
    std::vector<double> _N;
    std::vector<double> _Ngen;
    std::vector<std::vector<double> > _S;
    std::vector<std::vector<double> > _B;
    std::vector<double> _meanS;
    std::vector<double> _meanB;
    
    TRandom3 _random;
    int _nbins;
    int _index;
    bool _profile;
};

#endif

