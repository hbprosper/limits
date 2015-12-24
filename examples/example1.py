#!/usr/bin/env python
#-----------------------------------------------------------------------------
# example1.py
#
# Given
#      N     observed count
#      eff   signal efficiency * luminosity
#      deff  uncertainty
#      bkg   background estimate
#      dbkg  uncertainty in background
#
# compute a Bayesian upper limit on "cross section"
# Created June 2015 Les Houches 
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
def createInputs(filename, N, efl, defl, bkg, dbkg):
    # create input file for blimit.py
    out = open(filename, 'w')
    out.write('#----------------------------------------------------------\n')
    out.write('# Format:\n')
    out.write('#  bin1    bin2   ...\n')
    out.write('#  count1  count2 ...\n')
    out.write('#  efl1    efl2   ...\n')
    out.write('#  bkg1    bkg2   ...\n')
    out.write('#    :      :\n')
    out.write('# The first line is a header, which is followed by\n')
    out.write('# the counts, then the effective luminosities then\n')
    out.write('# the backgrounds\n')
    out.write('#\n')
    out.write('# In order to account for systematic uncertainty\n')
    out.write('# repeat the pair of effective luminosities and\n')
    out.write('# background lines with different random samplings.\n')
    out.write('#----------------------------------------------------------\n')
    out.write('%10s\n' % 'bin1')
    out.write('#counts\n')
    out.write('%10d\n' % int(N))
    if defl == 0 and dbkg == 0:
        out.write('#effective luminosities (efficiency x luminosity)\n')
        out.write('%10.4f\n' % efl)
        out.write('#backgrounds\n')
        out.write('%10.1f\n' % bkg)
    else:
        ntrial = 500
        gamma_efl, beta_efl = computeGammaConstants(efl, defl)
        gamma_bkg, beta_bkg = computeGammaConstants(bkg, dbkg)
        for ii in xrange(ntrial):
            out.write('#%d\n' % (ii+1))
            
            efl = gammavariate(gamma_efl, beta_efl)
            out.write('%10.4f\n' % efl)
            
            bkg = gammavariate(gamma_bkg, beta_bkg)
            out.write('%10.4f\n' % bkg)
    out.close()
#-----------------------------------------------------------------------------
def main():
    CL = 0.90
    # --------------------------------------
    # load limit codes
    # --------------------------------------
    gSystem.Load('liblimits')
    from ROOT import MultiPoisson, Bayes
         
    # --------------------------------------        
    # get results
    # --------------------------------------            
    N, efl, defl, bkg, dbkg = map(atof, sys.argv[1:])
    print "N   = %5.0f" % N
    print "eff = %10.4f, %-10.4f" % (efl, defl)
    print "bkg = %10.4f, %-10.4f" % (bkg, dbkg)
    print
    
    print "create inputs.dat"
    createInputs("inputs.dat", N, efl, defl, bkg, dbkg)
    
    # --------------------------------------        
    # create model
    # --------------------------------------
    print "create model"        
    model = MultiPoisson("inputs.dat")
    data  = model.counts()
    sigmamin =  0.0
    sigmamax = max(data)
    sigmamax += 5 * sqrt(sigmamax)

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
       python example1.py N efl defl bkg dbkg

       N     observed count
       efl   effective signal luminosity (signal efficiency X luminosity)
       defl  uncertainty in effective signal luminosity
       bkg   background estimate
       dbkg  uncertainty in background estimate
        ''')
    main()
    
except KeyboardInterrupt:
    print "ciao!"
    
