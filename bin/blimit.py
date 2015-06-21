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
    blimits.py input-file luminosity[=100/fb] CL[=0.90]

    input-file format:
    bin1    bin2   ...
    count1  count2 ...
    eff1    eff2   ...
    bkg1    bkg2   ...

    The first line is a header. Optionally, the last column is the
    luminosity. If the luminosity is given, we presume it is also
    sampled. Place the sampled luminosities at the end of each row
    of sampled efficiencies. Note: if the luminosity column is
    present, then the number of bins is one less than the number
    of columns.
    
    Each column (apart from the luminosity column) represents a bin,
    while each pair of lines of efficiencies and backgrounds contains
    random samplings of the predicted signal efficiencies and
    backgrounds. 

    This technique provides a simple, scalable, yet completely
    general, way to represent systematic uncertainty in predictions,
    without the need for assumptions about the manner in which the
    predictions are correlated across signals, backgrounds, and
    bins.
     
    The expected count (per bin) = xsec * eff * L + bkg

    Limits are set on the parameter "xsec". Note, however, that
    by setting L=1, and interpeting "eff" as the signal count,
    "xsec" can be interpreted as the signal strength "mu".
'''
CL = 0.90    # "confidence level"
#-----------------------------------------------------------------
def main():
    # --------------------------------------
    # load limit codes
    # --------------------------------------
    gSystem.Load('liblimits')
    from ROOT import MultiPoisson, Bayes

    argv = sys.argv[1:]
    filename = argv[0]
    if len(argv) < 2:
        luminosity = 100.0 # /fb
    else:
        luminosity = atof(argv[1])

    if len(argv) < 3:
        CL = 0.90
    else:
        CL = atof(argv[2])

    # --------------------------------------        
    # create model
    # --------------------------------------
    swatch = TStopwatch()
    swatch.Start()
    print "creating model"        
    model = MultiPoisson(filename, luminosity)
    print "\ttime =%8.3fs" % swatch.RealTime()            
    data     = model.counts()
    sigmamin =  0.0
    sigmamax = 20.0

    # --------------------------------------
    # compute limits
    # --------------------------------------
    swatch.Start()
    print "computing limit"        
    bayes = Bayes(model, data, sigmamin, sigmamax, CL)
    limit = bayes.quantile()
    print "=> limit: %6.1f fb (%2.0f%sCL)\n\ttime =%8.3fs" % \
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
