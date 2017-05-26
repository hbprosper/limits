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
    lowerlimit = wald.percentile(CLlow)
    upperlimit = wald.percentile(CLupp)
    print "=> central interval [%5.2f, %5.2f] (%4.1f%s) width = %5.2f" % \
      (lowerlimit, upperlimit, 100*CL, '%', upperlimit-lowerlimit)     

    CL = CLupper
    limit = wald.percentile(CL)
    print "=> upper limit: %5.2f (%2.0f%sCL)\t\ttime:  %8.3fs" % \
      (limit, 100*CL, '%', swatch.RealTime())

    # compute Z-value = sqrt(t0), where
    #              t0 = -2 ln [L(0)/L(poi_hat)]
    Z = wald.zvalue(0)
    print "=> Z-value:     %5.2f" % Z
    
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
    lowerlimit = bayes.percentile(CLlow)
    upperlimit = bayes.percentile(CLupp)
    print "=> central interval [%5.2f, %5.2f] (%4.1f%s) width = %5.2f" % \
      (lowerlimit, upperlimit, 100*CL, '%', upperlimit-lowerlimit)     

    CL = CLupper
    limit = bayes.percentile(CL)
    print "=> upper limit: %5.2f (%2.0f%sCL)\t\ttime:  %8.3fs" % \
      (limit, 100*CL, '%', swatch.RealTime())

    # compute Z-value, where
    #              Z = sign(lnB10)*sqrt(2*|lnB10|), B10 = p(s+b)/p(b)
    Z = bayes.zvalue(1)
    print "=> Z-value:     %5.2f" % Z

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
    
