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
	blimit.py inputs.dat        1.0 0.90
	number of bins 1
	counts
		0		1
	=> cross section range: [   0.0,   20.0]fb
	=> limit:    3.9 fb (90%CL)	time =   0.046s
	```

### Details

Given observed counts, signal efficiencies, and backgrounds, specified
in the file inputs.dat, blimit.py  computes Bayesian upper limits on
the signal cross section.

Usage:
```
    blimits.py input-file luminosity[=100/fb] CL[=0.90]
	```
The format of the input-file is:
	```
	bin1    bin2   ... 
    count1  count2 ...
    eff1    eff2   ...
    bkg1    bkg2   ...

    The first line is a header in which, optionally, the last column is the
    luminosity. If the luminosity is given, then a sampled value of
    the luminosity should be placed at the end of each row
    of sampled efficiencies. If the luminosity column is present, then
    the number of bins is presumed to be one fewer than the number of columns.
    
    Each column (apart from the luminosity column) represents a bin,
    while each pair of lines of efficiencies and backgrounds contains
    random samplings of the predicted signal efficiencies and
    backgrounds. The ensemble of efficiencies and backgrounds
    constitute the evidence-based prior with respect to which the
    likelihood is averaged.

    This technique provides a simple, scalable, yet completely
    general, way to represent systematic uncertainty in predictions,
    without the need for assumptions about the manner in which the
    predictions are correlated across signals, backgrounds, and
    bins.
     
    The expected count (per bin) = xsec * eff * L + bkg

    Limits are set on the parameter "xsec". Note, however, that
    by setting L=1, and interpeting "eff" as the signal count,
    "xsec" can be interpreted as the signal strength "mu".
	```
