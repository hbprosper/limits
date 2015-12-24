# limits
A simple standalone (multi-Poisson, multi-Poisson-gamma)  Bayes limit calculator

To build do
	make
  
To setup do
	source setup.sh
	or
	source setup.csh (for non-bash shells)

### Example 1
Usage:
```
python example1.py N eff deff bkg dbkg

     N         observed count
     eff       signal efficiency * luminosity estimate
     deff     uncertainty in signal efficiency estimate
     bkg     background estimate
     dbkg   uncertainty in background estimate
  ```
To test do
```
   cd examples
   python example1.py 1 1 0 0 0
   ```
  
Output:
```
   N     =     1
   eff   =     1.0000, 0.0000 
   bkg =     0.0000, 0.0000

   create inputs.dat
   create model
   compute limit
   => limit:    3.9 fb (90%CL)
      time:     0.059s
   ```


### Example 2
To test do
```
   cd examples
   python example2.py
   ```
This will run blimit.py (see below) on onebin.dat and threebin.dat.
Output:
```
	==> create model from onebin.dat <==
Wald
=> central interval [ 0.62,  1.22] (68.3%) width =  0.60
=> limit:  1.45 fb (95%CL)
   time:     0.024s

Bayes
=> central interval [ 0.64,  1.28] (68.3%) width =  0.65
=> limit:  1.52 fb (95%CL)
   time:     0.050s

	==> create model from threebin.dat <==
Wald
=> central interval [ 0.61,  1.20] (68.3%) width =  0.59
=> limit:  1.43 fb (95%CL)
   time:     0.024s

Bayes
=> central interval [ 0.65,  1.26] (68.3%) width =  0.61
=> limit:  1.50 fb (95%CL)
   time:     0.053s
```

### Details
Given observed counts, effective luminosities (efficiency*luminosity)
or predicted signal,
and backgrounds, specified in an file blimit.py  computes
Bayesian upper limits on the signal cross section (as well a
frequentist limit based on an asymptotic formula that makes use of the
Wald approximation (see "Asymptotic formulae for likelihood-based
tests of new physics", G. Cowan, K. Cranmer, E. Gross, and O. Vitells,
arXiv:1007.1727v3).

Usage:
```
    blimit.py input-file  [CL=0.95]
```
	
The format of the input-file is:
```
	bin1      bin2     ... 
    count1  count2 ...
    efl1       efl2      ...
    bkg1     bkg2    ...

    The first line is a header. Commented lines begin with a "#"
    
    Each column corresponds to a bin, while each pair of lines after
    the counts contains random samplings of the effective luminosities
	(or predicted signals) and backgrounds. 

    This technique provides a simple, scalable, yet completely
    general, way to represent systematic uncertainty in predictions,
    without the need to assume how the predictions are correlated
	across signals, backgrounds, and bins.
     
     The expected count (per bin) = sigma * efl + bkg
	 or mu * signal + bkg

	Limits are set on the parameter "sigma" or the signal strength "mu".
```
