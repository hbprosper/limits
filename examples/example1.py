#!/usr/bin/env python
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
def createBlimitInputs(N, eff, deff, bkg, dbkg):
    # create input file for blimit.py
    out = open('inputs.dat', 'w')
    out.write('#----------------------------------------------------------\n')
    out.write('# Format:\n')
    out.write('#  bin1    bin2   \n')
    out.write('#  count1  count2 ...\n')
    out.write('#  eff1    eff2   ...\n')
    out.write('#  bkg1    bkg2   ...\n')
    out.write('#    :      :\n')
    out.write('# The first line is a header, which is followed by\n')
    out.write('# the counts, which is followed by the efficiencies then\n')
    out.write('# the backgrounds\n')
    out.write('#\n')
    out.write('# In order to account for uncertainty in the efficiencies\n')
    out.write('# and the backgrounds (and optionally the luminosity, which\n')
    out.write('# can be the last column) repeat the pair of efficiency and\n')
    out.write('# background lines with different random samplings of the\n')
    out.write('# efficiencies and the backgrounds\n')
    out.write('#----------------------------------------------------------\n')
    out.write('%10s\n' % 'bin1')
    out.write('#counts\n')
    out.write('%10d\n' % int(N))
    if deff == 0 and dbkg == 0:
        out.write('#efficiencies\n')
        out.write('%10.4f\n' % eff)
        out.write('#backgrounds\n')
        out.write('%10.1f\n' % bkg)
    else:
        ntrial = 500
        gamma_eff, beta_eff = computeGammaConstants(eff, deff)
        gamma_bkg, beta_bkg = computeGammaConstants(bkg, dbkg)
        for ii in xrange(ntrial):
            out.write('#%d\n' % (ii+1))
            
            eff = gammavariate(gamma_eff, beta_eff)
            out.write('%10.4f\n' % eff)
            
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
    N, eff, deff, bkg, dbkg, luminosity = map(atof, sys.argv[1:])
    print "N   = %5.0f" % N
    print "eff = %10.4f, %-10.4f" % (eff, deff)
    print "bkg = %10.4f, %-10.4f" % (bkg, dbkg)
    print
    
    print "create inputs.dat"
    createBlimitInputs(N, eff, deff, bkg, dbkg)
    
    # --------------------------------------        
    # create model
    # --------------------------------------
    print "create model"        
    model = MultiPoisson("inputs.dat", luminosity)
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
    if len(argv) < 6:
        exit('''
    Usage:
       python example.py N eff deff bkg dbkg L

       N     observed count
       eff   signal efficiency estimate
       deff  uncertainty in signal efficiency estimate
       bkg   background estimate
       dbkg  uncertainty in background estimate
       L     integrated luminosity
        ''')
    main()
    
except KeyboardInterrupt:
    print "ciao!"
    
