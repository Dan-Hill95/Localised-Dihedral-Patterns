*These codes have either been written in 2021 by Dan J. Hill or adapted from codes by David J.B. Lloyd (University of Surrey) and Daniele Avitabile, "Numerical computation of coherent structures in spatially-extended neural networks", Second International Conference on Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.*

**For a more thorough tutorial in secant continuation and access to the original codes, see https://github.com/danieleavitabile/continuation-spatially-extended-systems-tutorial.**

Before accessing any of the available folders, you should initialise this code by running

Init.m

There are three types of folders:

**1. Passive folders - These are shared folders containing codes used by multiple routines. They are added to the path by running 'Init.m'.**

These include:

a) *Continuation_Codes* : containing the shared codes required to run path-following routines for numerical continuation.

b) *Equations* : containing the functions that define the discretised 2-3 Swift-Hohenberg equation, the bifurcation-tracking scheme, as well as the Galerkin truncated-Fourier system for approximating cellular patches.

**2. Standalone folders - These are able to be run, following the initialisation 'Init.m', without any prior data.**
			
These include:

a) *Radial_Solver* : Solves the 2-3 SH equation for localised radial patterns, tracks bifurcation curves, and computes linear stability.

b) *N4_Patches_N1Guess* : Solves the 2-3 SH equation for a localised D_{m} pattern with an explicit small truncation intial guess, tracks bifurcation curves, and computes linear stability.

c) *Patch_Match* : Solves the explicit algebraic matching for an approximate localised D_{m} pattern, and plots their planar profile for a given localisation.

**3. Piggyback folders - These are codes that require some prior data in order to be run. For all codes currently available, you should first run the codes in 'Radial_Solver', in order to obtain bifurcation data regarding localised radial patterns.**

These include:

a) *2_par_Bifurcation_Tracking* : Starting close to a change in linear stability, this code tracks the bifurcation point in parameter space, as a chosen parameter is varied.

b) *N1_Patches_Bif* : Starting close to a change in linear D_{m} stability, this code solves the coupled 2-3 SH equation for a D_{m}-perturbed localised radial pattern. These solutions form bifurcation curves that connect bifurcation points on the radial pattern's solution curve.

c) *N4_Patches_Bif* : Starting close to a change in linear stability, this code solves the 2-3 SH equation for a localised D_{m} pattern bifurcating from the localised radial solution.

For more information regarding the specific codes, see the ReadMe files in any of the 'Standalone' or 'Piggyback' folders.
