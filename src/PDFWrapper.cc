//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
// File: PDFWrapper.cc
// Description: Implements wrapper for RooAbsPdfs.
// 
// Created: June 11, 2010
//--------------------------------------------------------------
#include <iostream>
#include <cmath>
#include "TMath.h"
#ifdef __WITH_ROOFIT__
#include "RooDataSet.h"
#include "RooLinkedListIter.h"
#endif
#include "PDFWrapper.h"

ClassImp(PDFWrapper)

using namespace std;
//--------------------------------------------------------------
PDFWrapper::PDFWrapper()
  : PDFunction(),
    _data(vector<double>())
{}
#ifdef __WITH_ROOFIT__
PDFWrapper::PDFWrapper(RooAbsPdf& pdf, RooArgSet& obs, RooRealVar& poi) 
  : PDFunction(),
    _pdf(&pdf),
    _obs(&obs),
    _poi(&poi),
    _list(RooArgList(obs)),
    _data(vector<double>(_obs->getSize()))
{}

PDFWrapper::PDFWrapper(const PDFWrapper& other)
  : PDFunction(),
    _pdf(other._pdf),
    _obs(other._obs),
    _poi(other._poi),
    _list(other._list),
    _data(other._data)
{}

PDFWrapper::~PDFWrapper() 
{}

vector<double>&
PDFWrapper::generate(double poi)
{
  _poi->setVal(poi);
  const RooArgSet* row = _pdf->generate(*_obs, 1)->get();
  RooArgList d(*row);
  for(unsigned int c=0; c < _data.size(); c++)
    _data[c] = dynamic_cast<RooRealVar*>(&d[c])->getVal();
  delete row;
  return _data;
}

double 
PDFWrapper::operator() (std::vector<double>& data, double poi)
{
  if ( (int)data.size() == 0 ) return -1;
  if ( (int)data.size() != (int)_data.size() ) return -2;
  
  _poi->setVal(poi);
  for(unsigned int i=0; i < data.size(); i++)
    {
      _data[i] = data[i];
      RooRealVar* v = (RooRealVar*)(&_list[i]);
      v->setVal(_data[i]);
    }
  return _pdf->getVal();
}
#endif
