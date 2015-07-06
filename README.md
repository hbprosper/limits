# limits
A simple standalone (multi-Poisson)  Bayes limit calculator

To build do
	make
  
To setup do
	source setup.sh

### Example 1
Usage:
```
python example1.py N eff deff bkg dbkg L

     N     observed count
     eff   signal efficiency estimate
     deff  uncertainty in signal efficiency estimate
     bkg   background estimate
     dbkg  uncertainty in background estimate
     L     integrated luminosity
  ```
To test do
```
   cd example
   python example1.py 1 1 0 0 0 1
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

### Details
Given observed counts, effective luminosities (efficiency*luminosity),
and backgrounds, specified in the file inputs.dat, blimit.py  computes
Bayesian upper limits on the signal cross section (as well a
frequentist limit based on an asymptotic formula that makes use of the
Wald approximation (see "Asymptotic formulae for likelihood-based
tests of new physics", G. Cowan, K. Cranmer, E. Gross, and O. Vitells,
arXiv:1007.1727v3).

Usage:
```
    blimits.py input-file  CL[=0.90]
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
	and backgrounds. 

    This technique provides a simple, scalable, yet completely
    general, way to represent systematic uncertainty in predictions,
    without the need to assume how the predictions are correlated
	across signals, backgrounds, and bins.
     
    The expected count (per bin) = sigma * efl + bkg

    Limits are set on the parameter "sigma". Note, however, that
    by interpeting "efl" as the signal count, "sigma" can be
    interpreted as the signal strength "mu".
```
