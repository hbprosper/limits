#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File:        example3.py
# Description: Compute 
#              Example usage:
#                 example3.py input-file [xmin=0] [xmax=4] [CLupper=0.95]
#
# Created:     01-Jun-2016 HBP Bari
#              22-Jun-2016 HBP add ability to read histograms
#-----------------------------------------------------------------------------
import os, sys
from string import atof
from math import *
from ROOT import *
#-----------------------------------------------------------------------------
def main():
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
                
    print
    print "-"*50
    print "\t\texample3 - Wald"
    print "-"*50
    
    argv = sys.argv[1:]    
    filename = argv[0]
    if not os.path.exists(filename):
        sys.exit("** can't find file %s" % filename)

    if len(argv) < 2:
        mumin = 0.0
    else:
        mumin = atof(argv[1])

    if len(argv) < 3:
        mumax = 10.0
    else:
        mumax = atof(argv[2])
                
    if len(argv) < 4:
        CLupper = 0.95
    else:
        CLupper = atof(argv[3])
        
    # --------------------------------------        
    # create model
    # --------------------------------------
    print "=> inputs:  %s" % filename        
    model = MultiPoissonGamma(filename)    
    data  = model.counts()
    
    # --------------------------------------
    # compute Wald limits
    # --------------------------------------
    print "=> range:            [%8.2f, %8.2f]" % (mumin, mumax)
    
    wald = Wald(model, data, mumin, mumax)

    CL = 0.683
    CLlow = (1-CL)/2
    CLupp = (1+CL)/2
    lowerlimit = wald.percentile(CLlow)
    upperlimit = wald.percentile(CLupp)
    print "=> central interval: [%8.2f, %8.2f] (%4.1f%s CL)" % \
      (lowerlimit, upperlimit, 100*CL, '%')     

    CL = CLupper
    limit = wald.percentile(CL)
    print "=> upper limit:       %8.2f (%4.1f%s CL)" % (limit, 100*CL, '%')

    # compute Z-value = sqrt(t0), where
    #              t0 = 2 ln [L(poi_hat)/L(0)]
    # that is, twice the log of the likelihood of the best fit hypothesis
    # to the null hypothesis of no signal
    mu_hat = wald.estimate()
    dmu    = wald.uncertainty()
    
    Z = wald.zvalue(0)
    print "=> Z-value:           %8.2f" % Z

    print "=> mu_hat:            %8.2f +/- %-8.2f" % (mu_hat, dmu)
    
    relErr = mu_hat/dmu
    print "=> mu_hat/dmu (1):    %8.2f" % relErr
    
    relErr = (upperlimit+lowerlimit)/(upperlimit-lowerlimit)    
    print "=> mu_hat/dmu (2):    %8.2f" % relErr
#-----------------------------------------------------------------------------
try:
    argv = sys.argv[1:]
    if len(argv) < 1:
        exit('''
    Usage:
       example3.py filename [xmin=0] [xmax=10] [CL(upper)=0.95]
        ''')
    main() 
except KeyboardInterrupt:
    print "ciao!"
    
