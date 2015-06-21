#!/usr/bin/env python
#-----------------------------------------------------------------------------
import os, sys
from string import atof
from random import gammavariate
from math import *
from ROOT import gSystem
#-----------------------------------------------------------------------------
def computeGammaConstants(c, ec):
    # protect against negative counts
    # and zero uncertainty            
    if c  <= 0: c  = 1.e-6
    if ec <= 0: ec = 1.e-7
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
    
    # get results
    N, eff, deff, bkg, dbkg, L = map(atof, sys.argv[1:])
    print "N   = %5.0f" % N
    print "eff = %10.4f, %-10.4f" % (eff, deff)
    print "bkg = %10.4f, %-10.4f" % (bkg, dbkg)
    
    createBlimitInputs(N, eff, deff, bkg, dbkg)
    
    # compute limit
    cmd = 'blimit.py inputs.dat %10.1f %4.2f' % (L, CL)
    print cmd
    os.system(cmd)
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
    
