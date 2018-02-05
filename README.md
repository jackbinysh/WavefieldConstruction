A code which improves upon the "solid angle" method of vortex initialisaiton by consructing an entire wavefield.

Currently it can do the following:

Takes a file with a framed curve in it, and modifies the solid angle function to give a new phase field, phi, such that a chosen contour is tangent to the curve's framing in a tubular neighborhood of the curve.
Currently, *which* contour is hard coded in the function "ComputeSolidAngleFraming", under the name "phi0". It doesn't matter too much which it is, but phi0=0 might have trouble as its on the branch cut.

HOW TO USE ME:

name a .txt file in the "Constant.h" header, along with gridspacing, size etc. Make a .txt file of a curve with a framing - an example is given in the examples folder. It consists of point data, then, after a line with just an "A" on it, framing vectors for these points. Then run it - it will output "phi.vtk", the wavefield, and a curve, "Curvexxx.vtk", which is your input curve, with the target "A" framing and also the solid angle framing included.

compile as so:

export OMP_NUM_THREADS=16
g++ -fopenmp -O3 SolidAngleLink.cpp TriCubicInterpolator.cpp



