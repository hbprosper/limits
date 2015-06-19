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
    
    count1  count2 ....
    eff1    eff2 ...
    bkg1    bkg2...

    Each column represents a bin, while each pair of lines of
    efficiencies and backgrounds contains random samplings of
    the predicted signal efficiencies and backgrounds. This
    technique provides a simple, scalable, yet completely
    general, way to represent systematic uncertainty in predictions,
    without imposing assumptions about the manner in which the
    predictions are correlated across signals, backgrounds, and
    bins.
     
    expected count (per bin) = xsec * eff * L + bkg

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
    from ROOT import Bayes, CLsA, MultiPoisson

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
        
    # get data
    records = map(split, filter(lambda x: not (x[0] == '#' or x[0] == ''),
                                map(strip, open(filename).readlines())))

    # number of bins
    nbins = len(records[0])
    print "number of bins %d" % nbins
    
    # 1st row is row of observed counts    
    count = vector('double')(nbins)
    print "counts"
    total = 0
    for ii, N in enumerate(records[0]):
        count[ii] = atoi(N)
        print "\t%4d\t\t%d" % (ii, count[ii])
        total += count[ii]

    # --------------------------------------        
    # create model
    # --------------------------------------    
    model = MultiPoisson(count)        
    
    # get subsequent rows
    eff = vector('double')(nbins)
    bkg = vector('double')(nbins)
    records = records[1:]
    averageEff = 0.0
    averageBkg = 0.0
    nn = 0
    for ii in xrange(0, len(records), 2):
        nn += 1
        for jj in xrange(nbins):
            eff[jj] = atof(records[ii][jj])
            averageEff += eff[jj]
            
        for jj in xrange(nbins):
            bkg[jj] = atof(records[ii+1][jj])            
            averageBkg += bkg[jj]
            
        # add eff and bkg to model
        model.add(luminosity, eff, bkg)
        
    nn *= nbins
    averageEff /= nn
    averageBkg /= nn

    # --------------------------------------
    # set up cross section scan range
    # --------------------------------------
    a = averageEff * luminosity
    sigest = max(0.0, (total - averageBkg) / a)
    sigerr = sqrt(total)/a
    sigmamin = 0.0
    ii = int((sigest + 10 * sigerr)/10)
    sigmamax = (ii+1)*10
    print "=> cross section range: [%6.1f, %6.1f]fb\n" % \
      (sigmamin, sigmamax)

    # --------------------------------------
    # compute limits
    # --------------------------------------
    swatch = TStopwatch()
      
    # set up standalone Bayes limit calculator
    bayes = Bayes(model, count, sigmamin, sigmamax, CL)
    swatch.Start()    
    limitbayes = bayes.quantile()
    print "=> limit: %6.1f fb (%2.0f%sCL)\ttime =%8.3fs" % \
      (limitbayes, 100*CL, '%', swatch.RealTime())
#----------------------------------------------------------------------
try:
    argv = sys.argv[1:]
    if len(argv) < 1:
        exit(USAGE)
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
