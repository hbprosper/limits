# limits
A simple standalone (multi-Poisson)  Bayes limit calculator

To build do
	make
  
To setup do
	source setup.sh

### Example 1
	
Usage:
>	python example1.py N eff deff bkg dbkg L
>
>     N     observed count
>     eff   signal efficiency estimate
>     deff  uncertainty in signal efficiency estimate
>     bkg   background estimate
>     dbkg  uncertainty in background estimate
>     L     integrated luminosity

To test do
>   cd example
>   python example1.py 1 1 0 0 0 1
  
Output:
>
>	blimit.py inputs.dat        1.0 0.90
>	number of bins 1
>	counts
>		0		1
>	=> cross section range: [   0.0,   20.0]fb
>		
>	=> limit:    3.9 fb (90%CL)	time =   0.046s
>
