# 1d-hydrogel-modelling

This repository hosts the work made towards my summer research internship at the University of Birmingham studying 1d thermo-responsive hydrogels. Here you can find my final report (entitled <i>Shrinking and Swelling Dynamics in 1D Thermo-Responsive Hydrogels</i>) and a numerical scheme used to solve the non-linear PDE for porosity in a 1d gel.

There are 3 main .m files:

> `equilibriumCurve.m`: this produces a plot of $(T,\lambda)$ phase-space with the equilibrium
curve, spinodal and coexistence regions.

> `uniformTempDynamics.m`: this runs the solver for scenario where the temperature is abruptly
and uniformly changed throughout the gel. Try to keep temperatures in range 300-320 K, as 
otherwise `fsolve` often won't work. Need to pass the struct contained in `params.mat` as 
an argument to this function for it to work, together with the number of spatial steps, the 
time step, time range to solve for, the frequency of measurements of fluid fraction and gel length, 
whether you want a graphical output or not, and what file name to save under.

> `plotUniformTempDynamics.m`: this produces a 3d surface plot with time on the x-axis,
position on the y-axis and porosity on the z-axis.

All other files are utility functions and should be self-explanatory:

> `boundaryPhi.m` evaluates the porosity at the right-hand-boundary that satisfies the BC
there.

> `equilibriumT.m` calculates the equilibrium temperature for a given porosity.

> `firstDerivative.m` returns the first spatial derivative of an array, calculated using 
finite differences to a given accuracy - i'd suggest using accuracy=2.

> `params.mat` is the mat file which contains the parameters for the gel, e.g. the chi
function $\chi(T)$, $\Omega$ etc.


> `polymerVolume.m` is a file for evaluating the error in the calculated amount of polymer
in the gel at any given instant, used for validation.


