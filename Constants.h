//Constants.h
#include <string>

std::string knot_filename = "knotplot0_1001"; //assumed input filename format of "XXXXX.txt"
const int NumComponents = 1;   //No. components in the link

//OPTION - do you want the geometry of the input file to be exactly preserved, or can it be scaled to fit the box better
//1 to scale input file preserving the aspect ratio, 0 otherwise
#define PRESERVE_RATIOS 1  

// OPTION - what grid values do you want
//Grid points
const double h = 0.125;            //grid spacing
const int Nx = 200;   //No. points in x,y and z
const int Ny = 200;
const int Nz = 200;

// OPTION - how big should the knot be in the box, do you want it tilted or displaced?
//Size boundaries of knot (now autoscaled)
const double xmax = 8*Nx*h/10.0;
const double ymax = 8*Ny*h/10.0;
const double zmax = 8*Nz*h/10.0;

