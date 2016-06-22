#!/usr/bin/env python
#-----------------------------------------------------------------------------
# example2.py
#
# compute a Bayesian upper limit on signal strength.
#
# In this example, we assume a multi-Poisson-gamma model:
#     L(D|x, l, b) = prod_i=1^M PoissonGamma(N_i|x, l_i, dl_i, b_i, db_i)
# which is marginalized over a set of T points {(l_1,..l_M, b_1,..b_M)_i=j^T}
#
#-----------------------------------------------------------------------------
import os, sys
#-----------------------------------------------------------------------------
def main():
    os.system('blimit.py onebin.dat   0 4.0 0.95')
    os.system('blimit.py threebin.dat 0 4.0 0.95')     
#-----------------------------------------------------------------------------
try:
    main()    
except KeyboardInterrupt:
    print "ciao!"
    
