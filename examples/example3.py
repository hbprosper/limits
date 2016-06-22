#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File:        example3.py
# Description: Compute 
#              Example usage:
#                 example3.py input-file [xmin=0] [xmax=4] [CLupper=0.95]
#
# Created:     01-Jun-2016 HBP Bari
#              22-Jun-2016 HBP add ability to read histograms
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
    from ROOT import MultiPoissonGamma, Wald
        
    # --------------------------------------        
    # create model
    # --------------------------------------
    print "=> inputs:  %s" % filename        
    model = MultiPoissonGamma(filename)    
    data  = model.counts()
    
    # --------------------------------------
    # compute Wald limits
    # --------------------------------------
    print "=> range:            [%8.2f, %8.2f]" % (sigmamin, sigmamax)
    
    wald = Wald(model, data, sigmamin, sigmamax)

    CL = 0.683
    CLlow = (1-CL)/2
    CLupp = 1 - CLlow
    lowerlimit = wald.quantile(CLlow)
    upperlimit = wald.quantile(CLupp)
    print "=> central interval: [%8.2f, %8.2f] (%4.1f%s CL)" % \
      (lowerlimit, upperlimit, 100*CL, '%')     

    CL = CLupper
    limit = wald.quantile(CL)
    print "=> upper limit:       %8.2f (%4.1f%s CL)" % (limit, 100*CL, '%')

    # compute Z-value = sqrt(t0), where
    #              t0 = -2 ln [L(0)/L(poi_hat)]
    mu_hat = wald.estimate()
    dmu    = wald.uncertainty()
    
    Z = wald.zvalue(0)
    print "=> Z-value:           %8.2f" % Z

    print "=> mu_hat:            %8.2f +/- %-8.2f" % (mu_hat, dmu)
    
    relErr = mu_hat/dmu
    print "=> mu_hat/dmu (1):    %8.2f" % relErr
    
    relErr = (upperlimit+lowerlimit)/(upperlimit-lowerlimit)    
    print "=> mu_hat/dmu (2):    %8.2f" % relErr
#-----------------------------------------------------------------------------
try:
    argv = sys.argv[1:]
    if len(argv) < 1:
        exit('''
    Usage:
       example3.py filename [xmin=0] [xmax=4] [CL(upper)=0.95]
        ''')
    main() 
except KeyboardInterrupt:
    print "ciao!"
    
