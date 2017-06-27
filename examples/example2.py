#!/usr/bin/env python
#-----------------------------------------------------------------------------
# example2.py
#
# compute a Bayesian upper limit on signal strength.
#
# In this example, we assume a multi-Poisson-gamma model:
#     L(D|x, s, b) = prod_i=1^M PoissonGamma(N_i|x, s_i, ds_i, b_i, db_i)
# which is marginalized over a set of T points {(s_1,..s_M, b_1,..b_M)_i=j^T}
#
#-----------------------------------------------------------------------------
import os, sys
#-----------------------------------------------------------------------------
def main():
    print
    print "-"*50
    print "\t\texample2 - onebin, threebin"
    print "-"*50
            
    os.system('blimit.py onebin.dat   0 10.0 0.95')
    os.system('blimit.py threebin.dat 0 10.0 0.95')     
#-----------------------------------------------------------------------------
try:
    main()    
except KeyboardInterrupt:
    print "ciao!"
    
