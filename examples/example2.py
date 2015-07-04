#!/usr/bin/env python
#-----------------------------------------------------------------------------
import os, sys
from string import atof
from random import gammavariate
from math import *
from ROOT import gSystem, TFile, TStopwatch, kFALSE, kTRUE, vector
#-----------------------------------------------------------------------------
def computeGammaConstants(c, ec):
    # protect against negative counts
    # and zero uncertainty            
    if c  <= 0: c  = 1.e-3
    if ec <= 0: ec = 1.e-4
    k = c / ec
    k *= k
    gamma = (k+2 + sqrt((k+2)**2 - 4))/2
    beta  = (sqrt(c*c + 4*ec*ec) - c)/2
    return (gamma, beta)
#-----------------------------------------------------------------------------
def main():
    CL = 0.90
    # --------------------------------------
    # load limit codes
    # --------------------------------------
    gSystem.Load('liblimits')
    from ROOT import MultiPoissonGamma, Bayes
        
    # --------------------------------------        
    # get results
    # --------------------------------------            
    N, eff, deff, bkg, dbkg = map(atof, sys.argv[1:])[:5]
    print "N   = %5.0f" % N
    print "eff = %10.4f, %-10.4f" % (eff, deff)
    print "bkg = %10.4f, %-10.4f" % (bkg, dbkg)
    print
    
    # --------------------------------------        
    # create model
    # --------------------------------------
    print "create model"
    gamma_bkg, beta_bkg = computeGammaConstants(bkg, dbkg)    
    yy = vector('double')(1); yy[0] = gamma_bkg - 0.5
    b  = 1.0/beta_bkg
    
    gamma_eff, beta_eff = computeGammaConstants(eff, deff)
    xx = vector('double')(1); xx[0] = gamma_eff - 0.5
    a  = 1.0/beta_eff
    
    model = MultiPoissonGamma(yy, xx, b, a)
    data  = vector('double')(1, N)
    sigmamin =  0.0
    sigmamax = 20.0

    # --------------------------------------
    # compute limits
    # --------------------------------------
    swatch = TStopwatch()
    swatch.Start()
    
    print "compute limit"        
    bayes = Bayes(model, data, sigmamin, sigmamax, CL)
    limit = bayes.quantile()
    print "=> limit: %6.1f fb (%2.0f%sCL)\n   time:  %8.3fs" % \
      (limit, 100*CL, '%', swatch.RealTime())
#-----------------------------------------------------------------------------
try:
    argv = sys.argv[1:]
    if len(argv) < 5:
        exit('''
    Usage:
       python example2.py N eff deff bkg dbkg

       N     observed count
       eff   signal efficiency estimate
       deff  uncertainty in signal efficiency estimate
       bkg   background estimate
       dbkg  uncertainty in background estimate
        ''')
    main()
    
except KeyboardInterrupt:
    print "ciao!"
    
