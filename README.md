# limits
This is a simple (but reasonably general) standalone Bayesian limit calculator for analyses based on Poisson data, specifically, analyses that can be modeled using multi-Poisson or multi-Poisson-gamma likelihoods. For both models, the data comprise one or more observed counts. For the multi-Poisson model, we assume that a swarm of points (in 2 x _M_-dimensions, where _M_ is the number of counts) is available for the *parameters*, either the signals or the effective integrated luminosities (acceptance x efficiency x integrated luminosity) and backgrounds. The signal and background points are presumed to have been generated in another application through some sampling procedure. 

For the multi-Poisson-gamma model, in which the signal and background
parameters are marginalized over (that is, integrated out), we assume
that our knowledge of these parameters can be modeled using gamma
prior densities whose shape parameters are determined by identifying
each *estimate* of a parameter with the mode of the associated gamma
density and taking the *uncertainty* to be the gamma density's
standard deviation. Marginalization over systematic effects is
accomplished, as with the multi-Poisson model, by averaging over the
swarm of points in 2M-dimensions in which each point is the result of
a different random sampling of the systematic effects. Typically, the
sampling is done by identifying the sources of systematic error,
specifying algorithmically how each affects the objects whose
attributes ultimately feed into the signal and total background
bin-by-bin yields. 

For example, the true value of the *jet energy scale* (JES) is unknown; but we have an estimate of it and an associated uncertainty. Typically, both of these depend on the transverse momentum and rapidity of the jet. The probability density function (pdf) of the JES is often modeled as a Gaussian. For a given run of an analysis, a single random number _x_ is sampled from a zero mean unit width Gaussian and the transverse momentum of every jet in the sample of events is rescaled by 1 + _x_ * _sigma_, where _sigma_ is the standard deviation of the JES pdf for a given jet. The procedure is repeated a few hundred times, usually by running the multiple instances of the analysis in parallel. A similar procedure is followed for all other sources of systematic error. All such quantities should be sampled simultaneously. That way any complicated non-linear interactions induced in the signal and background yields by the systematic errors will be automatically taken into account without the need to assess how each of the  2 x _M_ yields are affected by these errors and what correlations, and what strength of correlation, exist between the yields due to the systematic errors. Furthermore, this procedure scales well with the number of systematic effects since all are sampled simultaneously.

## Setup
This package uses *ROOT* compiled with *mathmore*, which in turn requires the GNU Scientific Library (GSL), which is may be found at this website: http://www.gnu.org/software/gsl/

Once GSL is installed, install Root. You may have to tell ROOT the location of your GSL installation. To do so, just follow the README instructions that comes with ROOT. Then, to build _limits_ do
```
	make
```
  
To setup do
```
	source setup.sh
	or
	source setup.csh (for non-bash shells)
```
	
## Examples
	
### Example 1
Usage:
```
python example1.py N S dS B dB

	N	observed count
	S	signal (or effective luminosity = efficiency * luminosity)
	dS	uncertainty in signal
	B	background
	dB	uncertainty in background
  ```
To test do
```
   cd examples
   python example1.py 1 1 0 0 0
   ```
  
Output:
```
   N     =     1
   S     =     1.0000, 0.0000 
   B 	 =     0.0000, 0.0000

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
Given observed counts, predicted signals or effective integrated luminosities,
and backgrounds, specified in a text file, blimit.py  computes
Bayesian upper limits on the signal cross section (as well a
frequentist limit based on an asymptotic formula that makes use of the
Wald approximation (see "Asymptotic formulae for likelihood-based
tests of new physics", G. Cowan, K. Cranmer, E. Gross, and O. Vitells,
arXiv:1007.1727v3).

#### example1.py
This example uses the multi-Poisson model and a swarm of points in the
space of signals and backgrounds that
constitute a discrete representation of an evidence-based prior for
the latter. The multi-Poisson
likelihood is averaged over this prior. To use this model, you need to
have a model for how the signal and
background _estimates_ are probabilistically related to the
corresponding signal and background
_parameters_.  Consider, for example, the background model with parameter and
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
"s" is adequate, then just use example2.py.

The format of the input-file for the multi-Poisson model is:
```
number of bins
count1	count2 		...
number of sampled points
s1	s2         	...
b1	b2        	...
	
    Commented lines begin with a "#".
    
    Each column corresponds to a bin, while each pair of lines after
    the counts contains random samplings of the signal and background parameters. 

    This technique provides a simple, scalable, yet completely
    general, way to represent systematic uncertainty in predictions,
    without the need to make assumptions about how the predictions are correlated
    across signals, backgrounds, and bins.
     
     The expected count (per bin) is

	mu * s + b

	or

	sigma * l + b where l is the effective integrated luminosity.

	Limits are set on the parameter the signal strength
	"mu" or cross section "sigma" depending on what information is supplied. In the latter case, "l" is interpreted as the
	signal rather than the effective integrated luminosity.
```


#### example2.py
Unlike example1.py, which uses the
multi-Poisson model and an input file format in which the entries
after the counts are _parameter_ values, rather than estimates thereof,
*blimit.py* uses the _multi-Poisson-gamma_ model, which requires the
_estimates_ of the signals and backgrounds
together with their associated uncertainties. The idea here is to
simplify your life a little if the priors for the signal and backgrounds is of the form described above (see 
*example1.py*).

The input file for this example contain, not sampled parameters, but rather 
estimates of these quantities, along with their uncertainties, while
statistical dependencies  among
signals, backgrounds, and bins, is modeled by providing many sets of
estimates. The marginalization over these sets of estimates, which
models
integration over the systematic effects, is
approximated by Monte Carlo integration - basically, one
averages over the multi-Poisson-gamma model, while the marginalization
over the corresponding _parameters_ is done exactly. 

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
    blimit.py input-file  [xmin=0] [xmax=10] [CL=0.95]
```

The format of the input-file is:
```
number of bins
count1	count2 		...
number of samples
S1	S2         	...
dS1	dS2       	...
B1	B2        	...
dB1	dB2       	...
	: 
(See comments in example1.py above.)
```

## Using Histograms

When the number of bins is large (say > 10), it may be more convenient
to specify the counts and yields as ROOT histograms. The format is:
```
number of bins
data-root-filename  data-histogram-name
number of samples
signal-root-filename signal-histogram-name
background-root-filename background-histogram name
	:    :
```
Lines starting with # are treated as comments. The number of samples
is the number of pairs of signal/background files, which collectively
account for systematic uncertainties.
