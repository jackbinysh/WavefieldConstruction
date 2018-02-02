#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
using namespace std;

struct knotpoint
{
    double xcoord;   //position vector x coord
    double ycoord;   //position vector y coord
    double zcoord;   //position vector z coord
    double tx;       //tangent vector x component
    double ty;       //tangent vector y component
    double tz;       //tangent vector z component
    double ax;
    double ay;
    double az;
    double omegax;
    double omegay;
    double omegaz;
    double theta;
    double length;   //length of line
    double kappaNx;  //curvature vector x component
    double kappaNy;  //curvature vector x component
    double kappaNz;  //curvature vector x component
    double curvature;
    double twist;
};

struct knotcurve
{
    std::vector<knotpoint> knotcurve; //the actual data of the curve
    double length;   //total length of line
    double writhe;
    double twist;
};

struct Link
{
    std::vector<knotcurve> Components;
    int NumComponents;
    int NumPoints;
    double length;   //total length of line
    double writhe;
    double twist;
}; 

struct viewpoint
{
    double xcoord;
    double ycoord;
    double zcoord;
};

struct griddata 
{
    int Nx;
    int Ny;
    int Nz;
};

/*************************Functions for knot initialisation*****************************/

void InitialiseFromFile(struct Link& Curve);
void ScaleFunction(double *scale, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin);
void RefineCurve(struct Link& Curve);

/**********************Functions for curve geometry************************/

void ComputeLengths(struct Link& Curve);
void ComputeTangent(struct Link& Curve);
void ComputeKappaN(struct Link& Curve);
void ComputeWrithe(struct Link& Curve);
void ComputeTwist(struct Link& Curve);
void print_Curve(struct Link& Curve);
void print_B_phi(vector<double>&phi, const griddata &griddata);

/***********************Functions for outputting the solid angle*************************/

double SolidAngleCalc(const Link& Curve, const viewpoint& View);
void ComputeSolidAngle(vector<double>&phi, const struct Link& Curve);
void ComputeSolidAngleFraming(vector<double> &phi, Link &Curve);

/*************************General maths and integer functions*****************************/

// little inline guys
inline int incp(int i, int p, int N);    //increment i with p for periodic boundary
inline double x(int i,const griddata& griddata);
inline double y(int i,const griddata& griddata);
inline double z(int i,const griddata& griddata);
inline  int pt( int i,  int j,  int k,const griddata& griddata);      //convert i,j,k to single index
