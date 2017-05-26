#ifndef WALD_H
#define WALD_H
// ---------------------------------------------------------------------------
// File: Wald.h
// Description: compute limits using the statistic
//
//              q(poi) = -2*ln L(poi)/L(poi_hat)
//
//              and assuming the Wald approximation
//
//              q(poi) = (poi-poi_hat)^2 / var(poi_hat), poi_hat <= poi
//                    or 0 for poi_hat > poi.
//
//  See "Asymptotic formulae for likelihood-based tests of new physics",
//      G. Cowan, K. Cranmer, E. Gross, and O. Vitells, arXiv:1007.1727v3
// 
// Created: 04-Jun-2015 Harrison B. Prosper
//          26 May 2017 HBP include a base class (needed to allow
//                      polymorphism with ExpectedLimits class)
// ---------------------------------------------------------------------------
#include <vector>
#include "PDFunction.h"
#include "LimitCalculator.h"
// ---------------------------------------------------------------------------
/** Compute limits based on Wald approximation.
    <p>
    Approximate frequentist limits are computed based on the test statistic
    \f[
    q(\mu) = 2 * \ln L(\hat{\mu})/ L(\mu),
    \f]
    in which we are comparing the likelihood of the best fit hypothesis
    for \f$\mu\f$, namely \f$\hat{\mu}\f$, to the likelihood of alternative 
    hypothesis for the value of \f$\mu\f$.
    <p>
    This test statistic is useful, if the Wald approximation
    \f[
    q(\mu) \approx (\mu - \hat{\mu})^2 / \sigma^2, \\
    \quad\textrm{for } \hat{\mu} \leq \mu, \textrm{ and 0 otherwise}.
    \f]
    <p>
    For the special case of testing the null hypothesis \f$\mu = 0\f$ against
    the best fit hypothesis \f$\mu = \hat{\mu}\f$, the \f$Z-\textrm{value}\f$
    (the number of standard deviations the best fit is from the null) is 
    given \f$Z = \sqrt{q(0)}\f$.
    <p>
 See "Asymptotic formulae for likelihood-based tests of new physics",
      G. Cowan, K. Cranmer, E. Gross, and O. Vitells, arXiv:1007.1727v3
      for an instructive discussion.
 */
class Wald : public LimitCalculator
{
 public:
  ///
  Wald () {}

  /** Compute limits based on Wald approximation.
      @param model  - probability density function (pdf)
      @param data   - observed data
      @param poimin - minimum of parameter of interest
      @param poimax - maximum of parameter of interest
      @param CL     - confidence level
  */
  Wald(PDFunction& model,
       std::vector<double>& data,
       double poimin,
       double poimax,
       double CL=0.95);
  
  virtual ~Wald();

  PDFunction* pdf() {return _model;}
  
  /** Compute p-value given parameter of interest.
   */
  double operator()(double poi);
  
  /** Compute Z-value given parameter of interest using Z = sqrt[2*ln L(poi_hat)/L(0)].
   */
  double zvalue(double poi);
   /** Return best fit value.
    * This should be called after the interval calculation.
   */
  double estimate();
  /** Return uncertainty associated with estimate.
   * This should be called after the interval calculation.
   */
  double uncertainty();
  
  /// Compute percentile.
  double percentile(double CL=-1);

  ///
  void setRange(double poimin, double poimax) {_poimin=poimin;_poimax=poimax;}

  ///
  void setData(std::vector<double>& d);

  // For internal use.
  double fit(double guess=-1);
  double nll(double poi);

 private:
  PDFunction* _model; 
  std::vector<double> _data;

  double   _poimin;
  double   _poimax;
  double   _alpha;
  double   _poihat;
  double   _poierr;

  double   _f(double poi);
  
};


#endif
