#!/usr/bin/env python
#-----------------------------------------------------------------------------
import os, sys
from string import atof
from random import gammavariate
from math import *
from ROOT import gSystem, TFile, TStopwatch, kFALSE, kTRUE, vector
#-----------------------------------------------------------------------------
def createInputs(filename, N, efl, defl, bkg, dbkg):
    out = open(filename, 'w')
    out.write('#----------------------------------------------------------\n')
    out.write('# Format:\n')
    out.write('#  bin1    bin2   ...\n')
    out.write('#  count1  count2 ...\n')
    out.write('#  efl1    efl2   ...\n')
    out.write('#  defl1   defl2  ...\n')
    out.write('#  bkg1    bkg2   ...\n')
    out.write('#  dbkg1   dbkg2  ...\n')
    out.write('#    :      :\n')
    out.write('# The first line is a header, which is followed by\n')
    out.write('# the counts, then the effective luminosities & unc., then\n')
    out.write('# the backgrounds and uncertainties.\n')
    out.write('#\n')
    out.write('# In order to account for systematic uncertainty\n')
    out.write('# repeat the quadruplet of effective luminosities and\n')
    out.write('# background lines with different random samplings.\n')
    out.write('#----------------------------------------------------------\n')
    out.write('%10s\n' % 'bin1')
    out.write('#counts\n')
    out.write('%10d\n' % int(N))
    out.write('#effective luminosities (efficiency x luminosity)\n')
    out.write('%10.4f\n' % efl)
    out.write('%10.4f\n' % defl)
    out.write('#backgrounds\n')
    out.write('%10.4f\n' % bkg)
    out.write('%10.4f\n' % dbkg)
    out.close()
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
    N, efl, defl, bkg, dbkg = map(atof, sys.argv[1:])
    print "N   = %5.0f" % N
    print "eff = %10.4f, %-10.4f" % (efl, defl)
    print "bkg = %10.4f, %-10.4f" % (bkg, dbkg)

    print "create inputs2.dat"
    createInputs("inputs2.dat", N, efl, defl, bkg, dbkg)
    
    # --------------------------------------        
    # create model
    # --------------------------------------
    print "create model"        
    model = MultiPoissonGamma("inputs2.dat")    
    data  = model.counts()
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
       python example2.py N efl defl bkg dbkg
       
       N     observed count
       efl   effective signal luminosity (signal efficiency X luminosity)
       defl  uncertainty in effective signal luminosity
       bkg   background estimate
       dbkg  uncertainty in background estimate
        ''')
    main()
    
except KeyboardInterrupt:
    print "ciao!"
    
