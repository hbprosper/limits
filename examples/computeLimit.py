#!/usr/bin/env python
#-----------------------------------------------------------------------------
import os, sys
from string import atof
from random import gammavariate
from math import *
from ROOT import gSystem, TFile, TStopwatch, kFALSE, kTRUE, vector
#-----------------------------------------------------------------------------
def main():
    CL = 0.95
    filename = sys.argv[1]
    if not os.path.exists(filename):
        sys.exit("** can't find file %s" % filename)
        
    # --------------------------------------
    # load limit codes
    # --------------------------------------
    gSystem.Load('liblimits')
    from ROOT import MultiPoissonGamma, Bayes, Wald
        
    # --------------------------------------        
    # create model
    # --------------------------------------
    print "create model from %s" % filename        
    model = MultiPoissonGamma(filename)    
    data  = model.counts()
    sigmamin =  0.0
    sigmamax =  4.0

    # --------------------------------------
    # compute Wald limits
    # --------------------------------------
    swatch = TStopwatch()
    swatch.Start()

    print "\nWald"
    wald = Wald(model, data, sigmamin, sigmamax)

    CL = 0.683
    CLlow = (1-CL)/2
    CLupp = 1 - CLlow
    lowerlimit = wald.quantile(CLlow)
    upperlimit = wald.quantile(CLupp)
    print "=> central interval [%5.2f, %5.2f] (%4.1f%s) width = %5.2f" % \
      (lowerlimit, upperlimit, 100*CL, '%', upperlimit-lowerlimit)     

    CL = 0.95
    limit = wald.quantile(CL)
    print "=> limit: %5.2f fb (%2.0f%sCL)\n   time:  %8.3fs" % \
      (limit, 100*CL, '%', swatch.RealTime())


    # --------------------------------------
    # compute Bayes limits
    # --------------------------------------
    print "\nBayes"
    swatch = TStopwatch()
    swatch.Start()
    
    bayes = Bayes(model, data, sigmamin, sigmamax)
    
    CL = 0.683
    CLlow = (1-CL)/2
    CLupp = 1 - CLlow
    lowerlimit = bayes.quantile(CLlow)
    upperlimit = bayes.quantile(CLupp)
    print "=> central interval [%5.2f, %5.2f] (%4.1f%s) width = %5.2f" % \
      (lowerlimit, upperlimit, 100*CL, '%', upperlimit-lowerlimit)     

    CL = 0.95
    limit = bayes.quantile(CL)
    print "=> limit: %5.2f fb (%2.0f%sCL)\n   time:  %8.3fs" % \
      (limit, 100*CL, '%', swatch.RealTime())


       
#-----------------------------------------------------------------------------
try:
    argv = sys.argv[1:]
    if len(argv) < 1:
        exit('''
    Usage:
       python computeLimits.py filename
        ''')
    main()
    
except KeyboardInterrupt:
    print "ciao!"
    
