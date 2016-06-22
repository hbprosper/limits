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
// ---------------------------------------------------------------------------
#include <vector>
#include "PDFunction.h"
// ---------------------------------------------------------------------------
///
class Wald
{
 public:
  ///
  Wald () {}

  ///
  Wald(PDFunction& model,
       std::vector<double>& data,
       double poimin,
       double poimax,
       double CL=0.95);
  
  virtual ~Wald();
 
  /** Compute p-value given parameter of interest.
   */
  double operator()(double poi);

  /** Compute Z-value given parameter of interest using Z = sqrt[-2*ln L(0)/L(poi_hat)].
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
  
  /// Compute quantile.
  double quantile(double CL=-1);

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
  double   _Z;

  double   _f(double poi);
  
  ClassDef(Wald,1);
};


#endif
