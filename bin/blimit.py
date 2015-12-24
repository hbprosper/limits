#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File:        blimit.py
#
#              Example usage:
#                 blimit.py input-file [xmin=0] [xmax=4] [CLupper=0.95]
#
# Created:     17-Jun-2015 HBP Les Houches
#-----------------------------------------------------------------------------
import os, sys
from string import atof
from math import *
from ROOT import gSystem, TFile, TStopwatch, kFALSE, kTRUE, vector
#-----------------------------------------------------------------------------
def main():

    argv = sys.argv[1:]    
    filename = argv[0]
    if not os.path.exists(filename):
        sys.exit("** can't find file %s" % filename)

    if len(argv) < 2:
        sigmamin = 0.0
    else:
        sigmamin = atof(argv[1])

    if len(argv) < 3:
        sigmamax = 4.0
    else:
        sigmamax = atof(argv[2])
                
    if len(argv) < 4:
        CLupper = 0.95
    else:
        CLupper = atof(argv[3])
        
    # --------------------------------------
    # load limit codes
    # --------------------------------------
    gSystem.Load('liblimits')
    from ROOT import MultiPoissonGamma, Bayes, Wald
        
    # --------------------------------------        
    # create model
    # --------------------------------------
    print "\n\t==> create model from %s <==" % filename        
    model = MultiPoissonGamma(filename)    
    data  = model.counts()
    
    # --------------------------------------
    # compute Wald limits
    # --------------------------------------
    swatch = TStopwatch()
    swatch.Start()

    print "Wald\t\trange: [%8.1f,%8.1f]" % (sigmamin, sigmamax)
    
    wald = Wald(model, data, sigmamin, sigmamax)

    CL = 0.683
    CLlow = (1-CL)/2
    CLupp = 1 - CLlow
    lowerlimit = wald.quantile(CLlow)
    upperlimit = wald.quantile(CLupp)
    print "=> central interval [%5.2f, %5.2f] (%4.1f%s) width = %5.2f" % \
      (lowerlimit, upperlimit, 100*CL, '%', upperlimit-lowerlimit)     

    CL = CLupper
    limit = wald.quantile(CL)
    print "=> limit: %5.2f fb (%2.0f%sCL)\n   time:  %8.3fs" % \
      (limit, 100*CL, '%', swatch.RealTime())


    # --------------------------------------
    # compute Bayes limits
    # --------------------------------------
    print "\nBayes\t\trange: [%8.1f,%8.1f]" % (sigmamin, sigmamax)    
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

    CL = CLupper
    limit = bayes.quantile(CL)
    print "=> limit: %5.2f fb (%2.0f%sCL)\n   time:  %8.3fs" % \
      (limit, 100*CL, '%', swatch.RealTime())

#-----------------------------------------------------------------------------
try:
    argv = sys.argv[1:]
    if len(argv) < 1:
        exit('''
    Usage:
       blimit.py filename [CLupper=0.95]
        ''')
    main()
    
except KeyboardInterrupt:
    print "ciao!"
    
