#ifndef CLSA_H
#define CLSA_H
// ---------------------------------------------------------------------------
// File: CLsA.h
// Description: compute CLs asymptotic limit
// 
// Created: 07 Sep 2011 HBP
// Updated: 06 Mar 2014 HBP
//          13 Aug 2014 HBP - add possibility to cache Qsb and Qb
//          14 Aug 2014 HBP - use TH1 to smooth and interpolate
//          03 Jun 2015 HBP - divorce ClsA from CLs
// ---------------------------------------------------------------------------
#include <vector>
#include "PDFunction.h"
// ---------------------------------------------------------------------------
class CLsA
{
 public:
  CLsA () {}

  CLsA(PDFunction& model,
       std::vector<double>& data,
       double poimin,
       double poimax,
       double CL=0.95);
  
  virtual ~CLsA();
 
  /** Compute CLsA given parameter of interest.
   */
  double operator()(double poi);

  double CLsbA() {return _CLsbA;}
  
  double CLbA()  {return _CLbA;}
  
  /** Return cached observed q-statistic.
      This is calculated when CLs is calculated.
   */  
  double Qobs()  {return _Qobs;}

  double limit(double CL=-1);

  void setRange(double poimin, double poimax) {_poimin=poimin;_poimax=poimax;}
  void setData(std::vector<double>& d);

  // For internal use only
  double fit(double guess=-1);
  double nll(double poi);
  double Qstatistic(double poi);

 private:
  PDFunction* _model; 
  std::vector<double> _data;

  double   _poimin;
  double   _poimax;
  double   _alpha;

  double   _CLsbA;
  double   _CLbA;
  double   _Qobs;
  
  double   _poihat;
  double   _poierr;

  double   _f(double poi);
  double   _g(double poi);

  double   _Ddata;

  
  ClassDef(CLsA,1);
};


#endif
