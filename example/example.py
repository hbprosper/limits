#!/usr/bin/env python
#-------------------------------------------------------------------------------------
import os, sys
from string import atof
from random import gammavariate
from math import *
from ROOT import gSystem
#-------------------------------------------------------------------------------------
def computeGammaConstants(c, ec):
    # protect against negative counts            
    if c < 0: c = 0.0
    if ec > 0:
        k = (c / ec)**2
        gamma = (k + 2 + sqrt((k + 2)**2 - 4))/2
        beta  = (sqrt(c**2 + 4*ec**2) - c)/2
    else:
        # Note: the gamma reduces to an exponential function
        # exp(-x/beta) -> exp(c*x), with c = -1/beta
        beta  = 1.0e-4
        gamma = 1.0+beta
    return (gamma, beta)
#-------------------------------------------------------------------------------------
def createBlimitInputs(N, eff, deff, bkg, dbkg):
    # create input file for blimit.py
    out = open('inputs.dat', 'w')
    out.write('#-------------------------------------------------------------\n')
    out.write('# Format:\n')
    out.write('#  count1  count2 ...\n')
    out.write('#  eff1    eff2   ...\n')
    out.write('#  bkg1    bkg2   ...\n')
    out.write('#    :      :\n')
    out.write('# To account for uncertainty in the efficiency and background\n')
    out.write('# the efficiency and background lines are repeated with\n')
    out.write('# different samplings of the efficiency and background\n')
    out.write('#-------------------------------------------------------------\n')
    out.write('%10d\n' % int(N))

    if deff == 0 and dbkg == 0:
        out.write('%10.4f\n' % eff)
        out.write('%10.1f\n' % bkg)
    else:
        ntrial = 100
        alpha_eff, beta_eff = computeGammaConstants(eff, deff)
        alpha_bkg, beta_bkg = computeGammaConstants(bkg, dbkg)
        for ii in xrange(ntrial):
            eff = gammavariate(alpha_eff, beta_eff)
            out.write('%10.4f\n' % eff)
            bkg = gammavariate(alpha_bkg, beta_bkg)
            out.write('%10.1f\n' % bkg)
    out.close()
#-------------------------------------------------------------------------------------
def main():

    # get results
    N, eff, deff, bkg, dbkg, L = map(atof, sys.argv[1:])

    createBlimitInputs(N, eff, deff, bkg, dbkg)
    
    # compute limit
    cmd = 'blimit.py inputs.dat %10.2f' % L
    print cmd
    os.system(cmd)
#-------------------------------------------------------------------------------------
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
    
