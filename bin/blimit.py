#!/usr/bin/env python
#-----------------------------------------------------------------
# File:        blimit.py
#
#              Example usage:
#                 ./blimit.py input-file
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
    blimits.py input-file luminosity

    input-file format:
    count1  count2 ....
                          repeated
    eff1    eff2 ...
    bkg1    bkg2...
    
    expected count (per bin) = xsec * eff * L + bkg
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
    print "=> in cross section range: [%6.1f, %6.1f]fb\n" % \
      (sigmamin, sigmamax)

    # --------------------------------------
    # compute limits
    # --------------------------------------
    swatch = TStopwatch()
      
    # set up standalone Bayes limit calculator
    bayes = Bayes(model, count, sigmamin, sigmamax, CL)
    swatch.Start()    
    limitbayes = bayes.quantile()
    print "=> Bayes limit:               %6.1ffb\t%8.3fs" % \
      (limitbayes, swatch.RealTime())

    # set up standalone CLsA limit calculator
    CLSA  = CLsA(model, count, sigmamin, sigmamax, CL)
    swatch.Start()
    limitCLsA = CLSA.limit()
    print "=> CLs(A) limit:              %6.1ffb\t%8.3fs" % \
      (limitCLsA, swatch.RealTime())  
#----------------------------------------------------------------------
try:
    argv = sys.argv[1:]
    if len(argv) < 1:
        exit(USAGE)
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
