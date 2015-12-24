#!/usr/bin/env python
#-----------------------------------------------------------------------------
import os, sys
#-----------------------------------------------------------------------------
def main():
    os.system('blimit.py onebin.dat 0.95')
    os.system('blimit.py threebin.dat 0.95')     
#-----------------------------------------------------------------------------
try:
    main()    
except KeyboardInterrupt:
    print "ciao!"
    
