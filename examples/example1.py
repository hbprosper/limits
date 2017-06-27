#!/usr/bin/env python
#-----------------------------------------------------------------------------
# example1.py
#
# Given
#      N     observed count
#      s     signal
#            (or signal efficiency * acceptance * integrated luminosity)
#      ds    uncertainty in signal or effective integrated luminosity
#      b     background estimate
#      db    uncertainty in background
#
# compute a Bayesian upper limit on "cross section".
#
# In this example, we assume a multi-Poisson model:
#
#     L(D|x, l, b) = prod_i=1^M Poisson(N_i|x * l_i + b_i)
#
# which is marginalized over a set of T points {(l_1,..l_M, b_1,..b_M)_i=j^T}.
#
# This is more general than example 2, which assumes a particular model for the
# background prior, namely, a gamma density. 
#
# Created June 2015 Les Houches 
#-----------------------------------------------------------------------------
import os, sys
from ROOT import *
from math import *
from string import atof
from random import gammavariate
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
def createInputs(filename, N, S, dS, B, dB):
    # create input file for blimit.py
    out = open(filename, 'w')
    out.write('#----------------------------------------------------------\n')
    out.write('# Format:\n')
    out.write('#  M       number of bins\n')
    out.write('#  N1      N2 ...\n')
    out.write('#  K       number of sampled points\n')
    out.write('#  S1      S2   ...\n')
    out.write('#  B1      B2   ...\n')
    out.write('#    :      :\n')
    out.write('# The item is the number of bins. This is followed\n')
    out.write('# by the observed counts, the signal counts, then\n')
    out.write('# the background counts\n')
    out.write('#\n')
    out.write('# In order to account for uncertainties whatever their origin\n')
    out.write('# repeat the last pair of lines with different random samplings\n')
    out.write('# of signal and background counts.\n')
    out.write('#----------------------------------------------------------\n')
    out.write('# number of bins\n')
    out.write('\t1\n')
    out.write('# observed counts\n')
    out.write('\t%d\n' % int(N))
    if dS == 0 and dB == 0:
        out.write('# number of points\n')
        out.write('\t1\n')
        out.write('# signals or effective luminosities (eff x Lumi)\n')
        out.write('\t%10.4e\n' % S)
        out.write('# backgrounds\n')
        out.write('\t%10.4e\n' % B)
    else:
        ntrial = 400
        gamma_S, beta_S = computeGammaConstants(S, dS)
        gamma_B, beta_B = computeGammaConstants(B, dB)
        out.write('# number of sampled points\n')
        out.write('\t%d\n' % ntrial)
        for ii in xrange(ntrial):
            jj = ii+1
            
            if dS > 0:
                sig = gammavariate(gamma_S, beta_S)
            else:
                sig = 0.0
            out.write('# %4d signal\n' % jj)                
            out.write('\t%10.3f\n' % sig)
            
            if dB > 0:
                bkg = gammavariate(gamma_B, beta_B)
            else:
                bkg = 0.0
            out.write('# %4d background\n' % jj)                        
            out.write('\t%10.3f\n' % bkg)
    out.close()
#-----------------------------------------------------------------------------
def main():
    print
    print "-"*50
    print "\t\texample1 - single bin"
    print "-"*50    
    CL = 0.90
    # --------------------------------------
    # load limit codes
    # --------------------------------------
    if os.environ.has_key('LIMITS_PATH'):
        gSystem.AddDynamicPath("$LIMITS_PATH/lib")
        gSystem.Load('liblimits')
    else:
        sys.exit('''
    please do
        cd ..
        source setup.sh
    to define environment variable LIMITS_PATH
        ''')
        
    # --------------------------------------        
    # get results
    # --------------------------------------            
    N, S, dS, B, dB = map(atof, sys.argv[1:])
    print "N = %5.0f" % N
    print "S = %10.4f, %-10.4f" % (S, dS)
    print "B = %10.4f, %-10.4f" % (B, dB)
    print
    
    print "create inputs.dat"
    createInputs("inputs.dat", N, S, dS, B, dB)

    # --------------------------------------        
    # create model
    # --------------------------------------
    print "create model"        
    model = MultiPoisson("inputs.dat")
    data  = model.counts()
    mumin =  0.0
    mumax =  10.0

    # --------------------------------------
    # compute limits
    # --------------------------------------
    swatch = TStopwatch()
    swatch.Start()
    
    print "compute limit"        
    bayes = Bayes(model, data, mumin, mumax, CL)
    limit = bayes.percentile()
    print "=> limit: %6.1f fb @ %2.0f%s CL\n   time:  %8.3fs" % \
      (limit, 100*CL, '%', swatch.RealTime())
#-----------------------------------------------------------------------------
try:
    argv = sys.argv[1:]
    if len(argv) < 5:
        exit('''
    Usage:
       python example1.py N S dS B dB

       N     observed count
       S     signal estimate
       dS    uncertainty in signal estimate
       B     background estimate
       dB    uncertainty in background estimate

    example:
        python example1.py 0 1 0 0 0

        should give a 90%s upper limit of 2.3 for N=0
        ''' % '%')
    main()
    
except KeyboardInterrupt:
    print "ciao!"
    
