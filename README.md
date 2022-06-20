# IA2005
A simulation code of a multi-scale dynamic rupture propagation from Ide &amp; Aochi (JGR, 2005). The numerical scheme is based on 3D Boundary Integral Equation Method for a planar fault (Fukuyama &amp; Madariaga, 1995; 1998) and a renormalization technique (Aochi &amp; Ide, 2004). This provides a single rupture event, set by IA05.prm, among the 16384 possible earthquake scenarios. This code is suitable for understanding the theoretical concept of multi-scale earthquake rupture process by Ide &amp; Aochi (2005) or for utilizing the simulated models for further analyses (i.e. earthquake source statistics, earthquake scenarios for ground motion simulations). Do not hesitate to contact us if you have any idea.  

https://doi.org/10.5281/zenodo.1472238

## History
	26th October 2018: First updated.

## Programs provided here
	stoc1-fft-IA05-distrib.f90: main program
	ran1.f: subroutine for generating random number
	fourn-d.f: subroutine for FFT
	kernel31s_05Avril.f: subroutine for Green's function
	stoc1-fft-IA05-distrib.exe: executable on Windows, compiled with ifort on Windows10 for Win32 environement. 

## Compile
	> ifort stoc1-fft-IA05-distrib.f90 ran1.f fourn-d.f kernel31s_05Avril.f 

## Input
	IA05.prm: An integer in the first row, for the number of simulation, should be manually set from "case1-1.mag". Choose one number between (1, 16384).

## Notice for Utilisation
	See Ide and Aochi (2005) for the detail explanation of the given model parameters. 
	For scale interation (number=Y), grid size is ds=4x4^Y (m), where Y=(0,3); time step is dt=ds/(2xVp), where Vp=6 (km/s).
	Fault slip is implicitly asigned to the x-direction.
	The largest Y has an entire process of rupture propagation, while others are intermediate. 
	The case of Y=0 always exists, however the the initial crack is suddenly given at time=0 and it is visible.

## Outputs (Be careful, many files are generated on the current directory)
	XXXXX_stepYZZZ.dat (given simulation number XXXXX, iteration of scale Y(0 to 3), time step ZZZ): x-grid, y-grid, slip vel (m/s), slip (m), shear stress (MPa).
	outputXXXXXi.dat : indicating hypocenter choice.
	momentXXXXXY.dat: indicating seismic moment releae function at each iteration Y. (time (s), moment rate (N.m/s), moment (N.m), Mw)
	outputXXXXXf.dat : indicating seismic moment release function (time step, time (s), moment rate (N.m/s), moment (N.m)
	outputXXXXXY.dat : intermediate results at the end of each iteration Y.
	hetero.org: fault heterogeneous map.
	hoge2.dat: model paramters.

### References
	Ide, S. and H. Aochi (2005). Earthquakes as multiscale dynamic ruptures with heterogeneous fracture surface energy, J. Geophys. Res., 110, B11303, doi:10.1029/2004JB003591. (for theoretical framework. All the other references therein.)
	
	
### Applications (please give us your examples of applications).
	Aochi, H., J. Le Puth and S. Ide (2006). Attempts at using a dynamic rupture source model for ground-motion simulations, Eos Trans. AGU, 87(36), West. Pac. Geophys. Meet. Suppl., Abstract S11E-0160. (first example of ground motion simulation)
	Rohmer, J. and H. Aochi (2015). Impact of channel-like erosion patterns on the frequency-magnitude distribution of earthquakes, Geophys. J. Int., 202(1), 670-677, doi:10.1093/gji/ggv181
	Aochi, H. and A. Burnol (2018). Mechanism of the ML4.0 25th April 2016 earthquake in southwest of France in the vicinity of the Lacq gas field, J. Seismol., 22, 1139-1155, doi: 10.1007/s10950-018-9758-5. (for an example of application)
