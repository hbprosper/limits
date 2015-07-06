#!/usr/bin/env python
#-----------------------------------------------------------------
# File:        blimit.py
#
#              Example usage:
#                 ./blimit.py input-file luminosity[=100/fb] CL[=0.90]
#
# Created:     17-Jun-2015 HBP Les Houches
#-----------------------------------------------------------------
import os,sys,re
from string import *
from math import *
from time import sleep
from ROOT import gSystem, TFile, TStopwatch, kFALSE, kTRUE, vector
#-----------------------------------------------------------------
USAGE = '''
Usage:
    blimits.py input-file CL[=0.90]

    input-file format:
    bin1    bin2   ...
    count1  count2 ...
    efl1    efl2   ...
    bkg1    bkg2   ...

    The first line is a header. Commented lines begin with a "#"
    
    Each column corresponds to a given bin, while each pair of lines of
    effective luminosities and backgrounds contains random samplings of
    the predicted effective luminosities and backgrounds. 

    This technique provides a simple, scalable, yet completely
    general, way to represent systematic uncertainty in predictions,
    without the need for assumptions about the manner in which the
    predictions are correlated across signals, backgrounds, and
    bins.
     
    The expected count (per bin) = sigma * efl + bkg

    Limits are set on the parameter "sigma". Note, however, that
    by interpeting "efl" as the signal count, "sigma" can be
    interpreted as the signal strength "mu".

    blimits.py input-file CL[=0.90]
'''
CL = 0.90    # "confidence level"
#-----------------------------------------------------------------
def main():
    # --------------------------------------
    # load limit codes
    # --------------------------------------
    gSystem.Load('liblimits')
    from ROOT import MultiPoisson, Bayes, Wald

    argv = sys.argv[1:]
    filename = argv[0]
    if len(argv) < 2:
        CL = 0.90
    else:
        CL = atof(argv[1])

    # --------------------------------------        
    # create model
    # --------------------------------------
    swatch = TStopwatch()
    swatch.Start()
    print 
    print "create model"        
    model = MultiPoisson(filename)
    print "\ttime =%8.3fs" % swatch.RealTime()
    print           
    data     = model.counts()
    sigmamin =  0.0
    sigmamax = max(data)
    sigmamax += 5 * sqrt(sigmamax)

    # --------------------------------------
    # compute limits
    # --------------------------------------
    swatch.Start()
    print "computing Bayesian limit"        
    bayes = Bayes(model, data, sigmamin, sigmamax, CL)
    limit = bayes.quantile()
    print "=> limit: %6.1f fb (%2.0f%sCL)\n   time:  %8.3fs" % \
      (limit, 100*CL, '%', swatch.RealTime())
    print

    swatch.Start()
    print "computing limit (using an asymptotically valid formula given in"
    print '"Asymptotic formulae for likelihood-based tests of new physics"'
    print 'G. Cowan, K. Cranmer, E. Gross, and O. Vitells, arXiv:1007.1727v3)\n'
        
    wald = Wald(model, data, sigmamin, sigmamax, CL)
    limit = wald.quantile()
    print "=> limit: %6.1f fb (%2.0f%sCL)\n   time:  %8.3fs" % \
      (limit, 100*CL, '%', swatch.RealTime())      
#----------------------------------------------------------------------
try:
    argv = sys.argv[1:]
    if len(argv) < 1:
        exit(USAGE)
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
