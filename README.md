# limits
This is a simple (but still quite general) standalone Bayesian limit calculator for analyses based on Poisson data, specifically, analyses that can be modeled using multi-Poisson or multi-Poisson-gamma likelihoods. For both models, the data comprise one or more observed counts. For the multi-Poisson model, we assume that a swarm of points (in 2M-dimensions, where M is the number of counts) is available for the *parameters*, effective integrated luminosities (acceptance x efficiency x integrated luminosity) and backgrounds, which have been generated in another application through some sampling procedure. For the multi-Poisson-gamma model, in which the effective integrated luminosities and backgrounds are marginalized over, we assume that our knowledge of these parameters can be encoded in gamma prior densitues whose shape parameters are determined by identifying each *estimate* of a parameter with the mode of the associated gamma density and taking the *uncertainty* to be the gamma density's standard deviation. Marginalization over systematic effects is accomplished, as with the multi-Poisson model, by averaging over a swarm of points in 2M-dimensions in which each point corresponds to a different random sampling of the systematic effects.

### Setup
This package uses *Root* compiled with *mathmore*, which in turn requires the GNU Scientific Library (GSL), which is here http://www.gnu.org/software/gsl/

Once GSL is installed, install Root. You may have to tell Root the location of your GSL installation. To do so, just follow the README instructions that comes with Root. Then, to build _limits_ do
```
	make
```
  
To setup do
```
	source setup.sh
	or
	source setup.csh (for non-bash shells)
```
### Example 1
Usage:
```
python example1.py N l dl b db

     N         	observed count
     l       	signal efficiency * luminosity estimate
     dl     	uncertainty in signal efficiency estimate
     b		background estimate
     db   	uncertainty in background estimate
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
   bkg 	 =     0.0000, 0.0000

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
Wald		range: [     0.0,     4.0]
=> central interval [ 0.62,  1.22] (68.3%) width =  0.60
=> limit:  1.45 fb (95%CL)
   time:     0.024s

Bayes		range: [     0.0,     4.0]
=> central interval [ 0.66,  1.27] (68.3%) width =  0.61
=> limit:  1.51 fb (95%CL)
   time:     0.051s

	==> create model from threebin.dat <==
Wald		range: [     0.0,     4.0]
=> central interval [ 0.61,  1.20] (68.3%) width =  0.59
=> limit:  1.43 fb (95%CL)
   time:     0.024s

Bayes		range: [     0.0,     4.0]
=> central interval [ 0.65,  1.26] (68.3%) width =  0.60
=> limit:  1.49 fb (95%CL)
   time:     0.051s
```

### Details
Given observed counts, effective integrated luminosities (acceptable X
efficiency X integrated luminosity)
or predicted signals,
and backgrounds, specified in a text file, blimit.py  computes
Bayesian upper limits on the signal cross section (as well a
frequentist limit based on an asymptotic formula that makes use of the
Wald approximation (see "Asymptotic formulae for likelihood-based
tests of new physics", G. Cowan, K. Cranmer, E. Gross, and O. Vitells,
arXiv:1007.1727v3).

#### example1.py
This example uses the multi-Poisson model and a swarm of points in the
space of effective integrated luminosities and backgrounds that
constitute a discrete representation of an evidence-based prior for
the latter. The multi-Poisson
likelihood is averaged over this prior. To use this model, you need to
have a model for how the effective integrated luminosity and
background estimates are probabilistically related to the
corresponding effective integrated luminosity and background
parameters.  Consider, for example, the background model with parameter and
estimate "b" and "B", respectively. A reasonable model
for the background prior is 
```
	p(b | B ) = exp(-qb) (qb)^B / Gamma(B+1) 
```
where "q" is a scale factor that connects the background in a bin of a
 background control region to the corresponding bin in the signal
 region.  If "q" is not known precisely, it too will be modeled with
 some prior density.

In order to use example1.py, you need first to sample from p(b |
B) for every bin and every background because what goes into the input
file for the multi-Poisson model are not the estimates "B" themselves, but
rather their associated parameters "b". If the above model for "b" and
"l" is adequate, then just use example2.py.

The format of the input-file for the multi-Poisson model is:
```
	bin1      bin2     	...
    	count1  count2 		...
    	l1          l2         	...
    	b1         b2        	...
	
    The first line is a header. Commented lines begin with a "#".
    
    Each column corresponds to a bin, while each pair of lines after
    the counts contains random samplings of the (expected) effective luminosities
	(or predicted signals) and (expected) backgrounds. 

    This technique provides a simple, scalable, yet completely
    general, way to represent systematic uncertainty in predictions,
    without the need to make assumptions about how the predictions are correlated
	across signals, backgrounds, and bins.
     
     The expected count (per bin) is

	sigma * l + b

	or

	mu * l + b

	Limits are set on the parameter "sigma" or the signal strength
	"mu". In the latter case, "l" is interpreted as the expected
	signal rather than the expected integrated luminosity.
```


#### example2.py
Unlike example1.py, which uses the
multi-Poisson model and an input file format in which the entries
after the counts are _parameter_ values, rather than estimates thereof,
_blimit.py_ uses the _multi-Poisson-gamma_ model, which requires the
estimates of the effective integrated luminosities and backgrounds
together with their associated uncertainties. The idea here is to
simplify your life a little if the priors for the effective integrated
luminosities and backgrounds is of the form described above (see 
*example1.py*).

The input file for this example contain, not sampled parameters, but rather 
estimates of these quantities, along with their uncertainties, while
statistical dependencies  among
signals, backgrounds, and bins, is modeled by providing many sets of
estimates. The marginalization over these sets of estimates, which
models
integration of the systematic effects, is
approximated by Monte Carlo integration - basically, one
averages over the multi-Poisson-gamma model, while the marginalization
over over  corresponding parameters is done exactly. 

For example, you may have a histogram for the signal and one for the
background from which you extract estimates _x_ +/-_dx_ for each bin. These are
the values needed by _example2.py_. If you have many pairs of signal and
background histograms, each associated with a different random
sampling of systematic effects from, for example, the jet energy scale, jet energy
resolution, parton distribution functions, trigger efficiency, etc.,
this ensemble of histograms constitute the prior density that models
the uncertainty due to systematic effects.

Usage:
```
    blimit.py input-file  [xmin=0] [xmax=4] [CL=0.95]
```

The format of the input-file is:
```
	bin1      bin2     	...
    	count1  count2 		...
    	l1          l2         	...
	dl1        dl2       	...
    	b1         b2        	...
    	db1       db2       .	..
	
    The first line is a header. Commented lines begin with a "#"
    
    Each column corresponds to a bin, while each quadruplet of lines after
    the counts contains random samplings of the effective luminosities
	(or predicted signals) and associated uncertainties, followed by
	background estimates and uncertainties. 

    This technique provides a simple, scalable, yet completely
    general, way to represent systematic uncertainty in predictions,
    without the need to assume how the predictions are correlated
	across signals, backgrounds, and bins.
     
     The expected count (per bin) = sigma * l + b
	 or mu * l + b

	Limits are set on the parameter "sigma" or the signal strength
	"mu". (See comments in example1.py above.)
```
