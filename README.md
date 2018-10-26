# IA2005
A simulation code of a multi-scale dynamic rupture propagation from Ide &amp; Aochi (JGR, 2005) 

## History
	26th October 2018: First updated.

## Programs provided here
	stoc1-fft-IA05-distrib.f90: main program
	ran1.f: subroutine for generating random number
	fourn-d.f: subroutine for FFT
	kernel31s_05Avril.f: subroutine for Green's function

## Compile
	> ifort stoc1-fft-IA05-distrib.f90 ran1.f fourn-d.f kernel31s_05Avril.f 

## Input
	IA05.prm: An integer in the first row, for the number of simulation, should be manually set from "case1-1.mag".

## Output (Be careful, many files are generated on the current directory)
	XXXXX_stepYZZZ.dat (given simulation number XXXXX, iteration of scale Y(0 to 3), time step ZZZ): X grid, Y grid, slip vel (m/s), slip (m), shear stress (MPa).
	outputXXXXXi.dat : indicating hypocenter choice.
	momentXXXXXY.dat: indicating seismic moment releae function at each iteration Y. (time (s), moment rate (N.m/s), moment (N.m), Mw)
	outputXXXXXf.dat : indicating seismic moment release function (time step, time (s), moment rate (N.m/s), moment (N.m)
	outputXXXXXY.dat : intermediate results at the end of each iteration Y.
	hetero.org: fault heterogeneous map.
	hoge2.dat: model paramters.
