#include "SolidAngleLink.h"    
#include "Constants.h"    
#include "TriCubicInterpolator.h"
#include <math.h>
#include <string.h>

int main (void)
{
    Link Curve;
    vector<double> phi(Nx*Ny*Nz);
    griddata griddata;
    griddata.Nx = Nx;
    griddata.Ny = Ny;
    griddata.Nz = Nz;
    InitialiseFromFile(Curve);
    cout << "calculating the solid angle..." << endl;
    ComputeSolidAngle(phi,Curve);
    cout << "finding solid angle framing" << endl;
    ComputeSolidAngleFraming(phi,Curve);
    print_Curve(Curve);
    print_B_phi(phi,griddata);
    return 0;
}

// read in link from file
void InitialiseFromFile(Link& Curve)
{
    // first up, how many components are there?
    Curve.NumComponents = NumComponents;
    Curve.Components.resize(NumComponents);
    Curve.NumPoints = 0;

    /*  For recording max and min input values*/
    double maxxin = 0;
    double maxyin = 0;
    double maxzin = 0;
    double minxin = 0;
    double minyin = 0;
    double minzin = 0;

    for(int i=0; i<Curve.NumComponents; i++)
    {
        stringstream ss;
        string buff,filename,temp;

        ss.clear();
        ss.str("");
        if (Curve.NumComponents==1)
        {
            ss << knot_filename << ".txt";
        }
        else
        {
            ss << knot_filename << "_" << i <<  ".txt";
        }

        filename = ss.str();
        ifstream CurveInputStream;   //knot file(s)
        CurveInputStream.open(filename.c_str());
        // read in the points
        while(CurveInputStream.good())
        {
            double xcoord,ycoord,zcoord;    //temporary variables
            if(getline(CurveInputStream,buff))
            {
                ss.clear();
                ss.str("");
                ss << buff;
                if(buff.compare("A")==0) {break;}
                ss >> xcoord >> ycoord >> zcoord;
                knotpoint Point;
                Point.xcoord = xcoord;
                Point.ycoord = ycoord;
                Point.zcoord = zcoord;
                Curve.Components[i].knotcurve.push_back(Point);
                // track max and min input values
                if(xcoord>maxxin) maxxin = xcoord;
                if(ycoord>maxyin) maxyin = ycoord;
                if(zcoord>maxzin) maxzin = zcoord;
                if(xcoord<minxin) minxin = xcoord;
                if(ycoord<minyin) minyin = ycoord;
                if(zcoord<minzin) minzin = zcoord;
            }
            else break;
            // construct a point and put it on the curve
        }

        // now read in A
        int NP = Curve.Components[i].knotcurve.size();
        for(int s =0; s< NP; s++)
        {
            if(CurveInputStream.good())   //read in points for knot
            {
                double ax,ay,az;    //temporary variables
                if(getline(CurveInputStream,buff))
                {
                    ss.clear();
                    ss.str("");
                    ss << buff;
                    ss >> ax >> ay >> az;
                    Curve.Components[i].knotcurve[s].ax = ax;
                    Curve.Components[i].knotcurve[s].ay = ay;
                    Curve.Components[i].knotcurve[s].az = az;
                }
                else break;
                // construct a point and put it on the curve
            }
        }
        CurveInputStream.close();
        // keep track of how many total points are added to the link
        Curve.NumPoints += Curve.Components[i].knotcurve.size();
    }

    // now centre and scale to a standard size
    double midpoint[3];
    midpoint[0] = 0.5*(maxxin+minxin);
    midpoint[1] = 0.5*(maxyin+minyin);
    midpoint[2] = 0.5*(maxzin+minzin);
    double scale[3];
    ScaleFunction(scale,maxxin,minxin,maxyin,minyin,maxzin,minzin);
    for(int i=0; i<Curve.NumComponents; i++)
    {
        for(int s=0; s<Curve.Components[i].knotcurve.size(); s++)
        {
            Curve.Components[i].knotcurve[s].xcoord = scale[0]*(Curve.Components[i].knotcurve[s].xcoord - midpoint[0]);
            Curve.Components[i].knotcurve[s].ycoord = scale[1]*(Curve.Components[i].knotcurve[s].ycoord - midpoint[1]);
            Curve.Components[i].knotcurve[s].zcoord = scale[2]*(Curve.Components[i].knotcurve[s].zcoord - midpoint[2]);
        }
    }
    // basic geometry
    ComputeLengths(Curve);
    ComputeTangent(Curve);
    ComputeKappaN(Curve);
    // increase number of points on each link component -- the minimum number can be changed
    
    while (Curve.NumPoints < 600*Curve.NumComponents)
    {
        RefineCurve(Curve);
        ComputeLengths(Curve);
        ComputeTangent(Curve);
        ComputeKappaN(Curve);
        cout << "link has size " << Curve.NumPoints << endl;
    }
    
    ComputeWrithe(Curve);
    ComputeTwist(Curve);
}


void ComputeTwist(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for(int s=0; s<NP; s++)
        {
            // twist of this segment
            double ds = Curve.Components[i].knotcurve[s].length;
            double dxds = Curve.Components[i].knotcurve[s].tx;
            double dyds = Curve.Components[i].knotcurve[s].ty;
            double dzds = Curve.Components[i].knotcurve[s].tz;
            double bx = (Curve.Components[i].knotcurve[incp(s,1,NP)].ax - Curve.Components[i].knotcurve[s].ax)/ds;
            double by = (Curve.Components[i].knotcurve[incp(s,1,NP)].ay - Curve.Components[i].knotcurve[s].ay)/ds;
            double bz = (Curve.Components[i].knotcurve[incp(s,1,NP)].az - Curve.Components[i].knotcurve[s].az)/ds;
            Curve.Components[i].knotcurve[s].twist = (dxds*(Curve.Components[i].knotcurve[s].ay*bz - Curve.Components[i].knotcurve[s].az*by) + dyds*(Curve.Components[i].knotcurve[s].az*bx - Curve.Components[i].knotcurve[s].ax*bz) + dzds*(Curve.Components[i].knotcurve[s].ax*by - Curve.Components[i].knotcurve[s].ay*bx))/(2*M_PI*sqrt(dxds*dxds + dyds*dyds + dzds*dzds));
            //Add on writhe, twist
            Curve.Components[i].twist  += Curve.Components[i].knotcurve[s].twist*ds;
        }
    }
}

// scale the input link to a specified size
void ScaleFunction(double *scale, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin)
{
    int i;
    bool nonzeroheight[3];  //marker: true if this dimension has non zero height in stl file
    if(maxxin-minxin>0) { scale[0] = xmax/(maxxin-minxin); nonzeroheight[0] = true; }
    else { scale[0] = 1;  nonzeroheight[0] = false; }
    if(maxyin-minyin>0) { scale[1] = ymax/(maxyin-minyin); nonzeroheight[1] = true; }
    else { scale[1] = 1;  nonzeroheight[1] = false; }
    if(maxzin-minzin>0) { scale[2] = zmax/(maxzin-minzin); nonzeroheight[2] = true; }
    else { scale[2] = 1;  nonzeroheight[2] = false; }
#if PRESERVE_RATIOS
    double minscale=1000000000;
    int imin=3;
    for(i=0; i<3; i++)   //find minimum scale factor
    {
        if(scale[i] < minscale && nonzeroheight[i])
        {
            imin = i;
            minscale = scale[i];
        }
    }
    if(imin < 3)      //scale x,y, and z directions by same scale factor
    {
        for(i=0; i<3; i++) scale[i] = scale[imin];
    }
#endif
}

// computes the lengths of each segment of the link
void ComputeLengths(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for(int s=0; s<NP; s++)
        {
            double dx = (Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord - Curve.Components[i].knotcurve[s].xcoord);
            double dy = (Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord - Curve.Components[i].knotcurve[s].ycoord);
            double dz = (Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord - Curve.Components[i].knotcurve[s].zcoord);
            double deltas = sqrt(dx*dx+dy*dy+dz*dz);
            Curve.Components[i].knotcurve[s].length = deltas;
            Curve.Components[i].length += deltas;
        }
    }
}

// computes the unit tangent of the link
void ComputeTangent(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for(int s=0; s<NP; s++)        // central difference scheme
        {
            double dsp = Curve.Components[i].knotcurve[s].length;
            double dsm = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
            Curve.Components[i].knotcurve[s].tx = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].xcoord - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].xcoord;
            Curve.Components[i].knotcurve[s].ty = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].ycoord - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].ycoord;
            Curve.Components[i].knotcurve[s].tz = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].zcoord - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].zcoord;
        }
    }
}

// computes the derivative of the unit tangent -- kappa N
void ComputeKappaN(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for(int s=0; s<NP; s++)    // central difference scheme
        {
            double dsp = Curve.Components[i].knotcurve[s].length;
            double dsm = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
            double kappaNx = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].tx + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].tx - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].tx;
            double kappaNy = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].ty + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].ty - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].ty;
            double kappaNz = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].tz + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].tz - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].tz;
            Curve.Components[i].knotcurve[s].kappaNx = kappaNx;
            Curve.Components[i].knotcurve[s].kappaNy = kappaNy;
            Curve.Components[i].knotcurve[s].kappaNz = kappaNz;
            // no longer need this -- could remove
            Curve.Components[i].knotcurve[s].curvature = sqrt(kappaNx*kappaNx + kappaNy*kappaNy + kappaNz*kappaNz);
        }
    }
}

// double the number of points on the curve -- simple interpolation
void RefineCurve(Link& Curve)
{
    Link NewCurve;
    NewCurve.NumComponents = NumComponents;
    NewCurve.Components.resize(NumComponents);
    NewCurve.NumPoints = 0;

    for (int i=0; i<Curve.NumComponents; i++) // run over the components
    {
        int NP = Curve.Components[i].knotcurve.size();
        for (int s=0; s<NP; s++) // run over the points of each component
        {
            knotpoint Point;
            // keep old point
            Point.xcoord = Curve.Components[i].knotcurve[s].xcoord;
            Point.ycoord = Curve.Components[i].knotcurve[s].ycoord;
            Point.zcoord = Curve.Components[i].knotcurve[s].zcoord;
            NewCurve.Components[i].knotcurve.push_back(Point);
            // create new point
            double ds = 0.5*Curve.Components[i].knotcurve[s].length;
            double x1 = Curve.Components[i].knotcurve[s].xcoord + ds*Curve.Components[i].knotcurve[s].tx + 0.5*ds*ds*Curve.Components[i].knotcurve[s].kappaNx;
            double x2 = Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord - ds*Curve.Components[i].knotcurve[incp(s,1,NP)].tx + 0.5*ds*ds*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNx;
            Point.xcoord = 0.5*(x1+x2);
            double y1 = Curve.Components[i].knotcurve[s].ycoord + ds*Curve.Components[i].knotcurve[s].ty + 0.5*ds*ds*Curve.Components[i].knotcurve[s].kappaNy;
            double y2 = Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord - ds*Curve.Components[i].knotcurve[incp(s,1,NP)].ty + 0.5*ds*ds*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNy;
            Point.ycoord = 0.5*(y1+y2);
            double z1 = Curve.Components[i].knotcurve[s].zcoord + ds*Curve.Components[i].knotcurve[s].tz + 0.5*ds*ds*Curve.Components[i].knotcurve[s].kappaNz;
            double z2 = Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord - ds*Curve.Components[i].knotcurve[incp(s,1,NP)].tz + 0.5*ds*ds*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNz;
            Point.zcoord = 0.5*(z1+z2);

            // interpolate A values for this new point - super quick and dirty, just take the average:
            Point.ax = 0.5*(Curve.Components[i].knotcurve[s].ax + Curve.Components[i].knotcurve[incp(s,1,NP)].ax);
            Point.ay = 0.5*(Curve.Components[i].knotcurve[s].ay + Curve.Components[i].knotcurve[incp(s,1,NP)].ay);
            Point.az = 0.5*(Curve.Components[i].knotcurve[s].az + Curve.Components[i].knotcurve[incp(s,1,NP)].az);
            NewCurve.Components[i].knotcurve.push_back(Point);

        }
        NewCurve.NumPoints += NewCurve.Components[i].knotcurve.size();
    }
    Curve = NewCurve;


}

// computes the writhe of each link component 
void ComputeWrithe(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        double Wr = 0.0;
        int NP = Curve.Components[i].knotcurve.size();

        for (int s=0; s<NP; s++) // running over the knot
        {
            double a1x = Curve.Components[i].knotcurve[s].xcoord;
            double a1y = Curve.Components[i].knotcurve[s].ycoord;
            double a1z = Curve.Components[i].knotcurve[s].zcoord;
            double t1x = Curve.Components[i].knotcurve[s].tx;
            double t1y = Curve.Components[i].knotcurve[s].ty;
            double t1z = Curve.Components[i].knotcurve[s].tz;
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);

            for (int t=s+1; t<NP; t++) // run over all points ahead of s
            {
                double a2x = Curve.Components[i].knotcurve[t].xcoord;
                double a2y = Curve.Components[i].knotcurve[t].ycoord;
                double a2z = Curve.Components[i].knotcurve[t].zcoord;
                double t2x = Curve.Components[i].knotcurve[t].tx;
                double t2y = Curve.Components[i].knotcurve[t].ty;
                double t2z = Curve.Components[i].knotcurve[t].tz;
                double dt = 0.5*(Curve.Components[i].knotcurve[t].length+Curve.Components[i].knotcurve[incp(t,-1,NP)].length);

                double dist = sqrt((a1x-a2x)*(a1x-a2x)+(a1y-a2y)*(a1y-a2y)+(a1z-a2z)*(a1z-a2z));
                Wr += ds*dt*((a1x-a2x)*(t1y*t2z-t1z*t2y)+(a1y-a2y)*(t1z*t2x-t1x*t2z)+(a1z-a2z)*(t1x*t2y-t1y*t2x))/(dist*dist*dist);
            }
        }
        Wr /= 2.0*M_PI;
        Curve.Components[i].writhe=Wr;
    }
}

// computes the solid angle the link presents from a viewpoint in its complement
// this version sweeps the solid angle out by moving the link in from being initially asymptotically far away
double SolidAngleCalc(const Link& Curve, const viewpoint& View)
{
    double totalomega = 0;
    double mindist = 1000.0; // scroll waves
    for(int i=0; i<Curve.NumComponents; i++)
    {
        double Integral = 0;
        int NP = Curve.Components[i].knotcurve.size();

        // first define and check a choice of asymptotic direction
        double ndotnmin = 1.0;
        double ndotnmax = -1.0;
        int smin;
        // define the asymptotic direction -- ninfty -- z-axis by default (in lower half space)
        double ninftyx = 0.0;
        double ninftyy = 0.0;
        double ninftyz = 1.0;
        if (View.zcoord>0) {ninftyz = -1.0;} // minus z in the upper half space
        for (int s=0; s<NP; s++)
        {
            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
            double ndotninfty = viewz*ninftyz/dist;
            if (ndotninfty<ndotnmin) {ndotnmin = ndotninfty; smin = s;}
            if (ndotninfty>ndotnmax) {ndotnmax = ndotninfty;}
        }
        if (ndotnmin < -0.98) // check if a threshold is exceeded -- value can be changed
        {
            if (ndotnmax < 0.98) {ninftyz = -ninftyz;} // flip direction
            else                                       // unless another threshold is exceeded -- value can be changed
            {
                ninftyz = 0.0;
                ninftyx = Curve.Components[i].knotcurve[smin].ty;    // set an orthogonal direction -- not guaranteed to be a good choice
                ninftyy = -Curve.Components[i].knotcurve[smin].tx;
                double norm = sqrt(ninftyx*ninftyx + ninftyy*ninftyy);
                ninftyx /= norm;
                ninftyy /= norm;
                //	  if (View.xcoord>0) {ninftyx = -ninftyx;} // could be picky about signs -- old code here, beware !!
            }
        }

        // scroll wave
        double mindistComp = 1000.0;
        // now for the calculation
        for (int s=0; s<NP; s++)
        {
            // define the view vector -- n = (Curve - View)/|Curve - View|
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
            double ndotninfty = viewx*ninftyx + viewy*ninftyy + viewz*ninftyz;
            double tx = Curve.Components[i].knotcurve[s].tx;
            double ty = Curve.Components[i].knotcurve[s].ty;
            double tz = Curve.Components[i].knotcurve[s].tz;
            // trapezium rule quadrature
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
            // and here's the integrand
            Integral += (ds/dist)*(ninftyz*(ty*viewx-tx*viewy)+ninftyx*(tz*viewy-ty*viewz)+ninftyy*(tx*viewz-tz*viewx))/(dist + ndotninfty);
            // scroll wave
           // if (dist < mindistComp) {mindistComp = dist;}
        }

        totalomega += Integral;
        // scroll wave
       // if (mindistComp < mindist) {mindist = mindistComp;}
    }
    //scroll wave
    // put in the interval [0,4pi]
    while(totalomega>4*M_PI) totalomega -= 4*M_PI;
    while(totalomega<0) totalomega += 4*M_PI;

    //  return totalomega;
   // return 0.5*totalomega-(16.0*M_PI/Nx)*mindist; // scroll wave
    return totalomega; 
}

void ComputeSolidAngle(vector<double>& phi, const Link& Curve)
{
    griddata griddata;
    griddata.Nx = Nx;
    griddata.Ny = Ny;
    griddata.Nz = Nz;
    double SolidAngle;
    viewpoint Point;

    for(int k=0; k<Nz; k++)
    {
        for(int j=0; j<Ny; j++)
        {
            for(int i=0; i<Nx; i++)
            {
                int n = pt(i,j,k,griddata);
                Point.xcoord = x(i,griddata);
                Point.ycoord = y(j,griddata) ;
                Point.zcoord = z(k,griddata);

                SolidAngle = SolidAngleCalc(Curve,Point);
                // put in the interval [0,4pi]
                while(SolidAngle>4*M_PI) SolidAngle -= 4*M_PI;
                while(SolidAngle<0) SolidAngle += 4*M_PI;
                phi[n]= SolidAngle/2;

            }
        }
    }
}

// find the difference between the framing we have given the curve, and the framing provided by the solid angle function.
void ComputeSolidAngleFraming(vector<double>&phi,Link&  Curve)
{
    // Im going to want phi values around the filament. for this, we construct an interpolator object
    likely::TriCubicInterpolator interpolatedphi(phi,h,Nx,Ny,Nz);

    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();

        for (int s=0; s<NP; s++)
        {
            // grab the point, compute the framing
            double xcoord= Curve.Components[i].knotcurve[s].xcoord;
            double ycoord= Curve.Components[i].knotcurve[s].ycoord;
            double zcoord= Curve.Components[i].knotcurve[s].zcoord;
            double ax= Curve.Components[i].knotcurve[s].ax;
            double ay= Curve.Components[i].knotcurve[s].ay;
            double az= Curve.Components[i].knotcurve[s].az;
            double tx= Curve.Components[i].knotcurve[s].tx;
            double ty= Curve.Components[i].knotcurve[s].ty;
            double tz= Curve.Components[i].knotcurve[s].tz;
            double tcax = ty*az-tz*ay;
            double tcay = tz*ax-tx*az;
            double tcaz = tx*ay-ty*ax;

            // okay, we have the frame at this point. lets walk in a small circle
            // around the filament and just find the rough value of phi closest to phi0.
            double phi0 = 3.14;
            int NumTestPoints=100;
            // impossible values for phi and theta
            double mintheta = -10;
            double minphi = 100;
            double finalvx=0 ;
            double finalvy=0;
            double finalvz=0 ;
            double r = 1;
            for(int q=0;q<NumTestPoints;q++)
            {
                double theta = ( ((double)q) / (double(NumTestPoints)) )*2*M_PI;
                double vx = r*(cos(theta)*ax + sin(theta)*tcax);
                double vy = r*(cos(theta)*ay + sin(theta)*tcay);
                double vz = r*(cos(theta)*az + sin(theta)*tcaz);

                double testphi = interpolatedphi(xcoord+vx,ycoord+vy,zcoord+vz);
                if(fabs(testphi-phi0) < fabs(minphi-phi0))
                {
                    minphi = testphi;
                    mintheta = theta;
                    finalvx =vx;
                    finalvy =vy;
                    finalvz =vz;
                }
            }
            Curve.Components[i].knotcurve[s].theta = mintheta;
            Curve.Components[i].knotcurve[s].omegax = finalvx;
            Curve.Components[i].knotcurve[s].omegay = finalvy;
            Curve.Components[i].knotcurve[s].omegaz = finalvz;
        }

    }
}

void print_B_phi( vector<double>&phi, const griddata& griddata)
{
    int Nx = griddata.Nx;
    int Ny = griddata.Ny;
    int Nz = griddata.Nz;
    int i,j,k,n;
    string fn = "phi.vtk";

    ofstream Bout (fn.c_str());

    Bout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET STRUCTURED_POINTS\n";
    Bout << "DIMENSIONS " << Nx << ' ' << Ny << ' ' << Nz << '\n';
    Bout << "ORIGIN " << x(0,griddata) << ' ' << y(0,griddata) << ' ' << z(0,griddata) << '\n';
    Bout << "SPACING " << h << ' ' << h << ' ' << h << '\n';
    Bout << "POINT_DATA " << Nx*Ny*Nz << '\n';
    Bout << "SCALARS Phi float\nLOOKUP_TABLE default\n";
    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Ny; j++)
        {
            for(i=0; i<Nx; i++)
            {
                n = pt(i,j,k,griddata);
                Bout << phi[n] << '\n';
            }
        }
    }
    Bout.close();
}


void print_Curve(struct Link& Curve)
{
    for( int c=0; c < (Curve.Components.size()) ; c++)
    {

        /***Write values to file*******/
        stringstream ss;

        ss.str("");
        ss.clear();

        ss << "Curve" << c << ".vtk";
        ofstream knotout (ss.str().c_str());

        int i;
        int n = Curve.Components[c].knotcurve.size();

        knotout << "# vtk DataFile Version 3.0\nKnot\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        knotout << "POINTS " << n << " float\n";

        for(i=0; i<n; i++)
        {
            knotout << Curve.Components[c].knotcurve[i].xcoord << ' ' << Curve.Components[c].knotcurve[i].ycoord << ' ' << Curve.Components[c].knotcurve[i].zcoord << '\n';
        }

        knotout << "\n\nCELLS " << n << ' ' << 3*n << '\n';

        for(i=0; i<n; i++)
        {
            knotout << 2 << ' ' << i << ' ' << incp(i,1,n) << '\n';
        }

        knotout << "\n\nCELL_TYPES " << n << '\n';

        for(i=0; i<n; i++)
        {
            knotout << "3\n";
        }

        knotout << "\n\nPOINT_DATA " << n << "\n\n";

        knotout << "\nVECTORS A float\n";
        for(i=0; i<n; i++)
        {
            knotout << Curve.Components[c].knotcurve[i].ax << ' ' << Curve.Components[c].knotcurve[i].ay << ' ' << Curve.Components[c].knotcurve[i].az << '\n';
        }

        knotout << "\nVECTORS Omega float\n";
        for(i=0; i<n; i++)
        {
            knotout << Curve.Components[c].knotcurve[i].omegax << ' ' << Curve.Components[c].knotcurve[i].omegay << ' ' << Curve.Components[c].knotcurve[i].omegaz << '\n';
        }
        knotout.close();
    }
}


inline int incp(int i, int p, int N)    //increment i with p for periodic boundary
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}
inline double x(int i,const griddata& griddata)
{
    return (i+0.5-griddata.Nx/2.0)*h;
}
inline double y(int i,const griddata& griddata)
{
    return (i+0.5-griddata.Ny/2.0)*h;
}
inline double z(int i,const griddata& griddata)
{
    return (i+0.5-griddata.Nz/2.0)*h;
}
inline  int pt( int i,  int j,  int k,const griddata& griddata)       //convert i,j,k to single index
{
    return (i*griddata.Ny*griddata.Nz+j*griddata.Nz+k);
}
