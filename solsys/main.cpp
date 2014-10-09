#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <planet.h>

using namespace std;


const double pi= 3.1415926;
const double G=4*pi*pi;
double forcex(double x,double y);
double forcey(double x,double y);
void write(double *z, double *y, int n, char *file);

//time - year, mass - sun mass, distance - AU

void RK4(double *vx,double *vy,double *x,double *y,int n,double h);

int main()
{
   double m_S=1;
   double *vx;
   double *vy;
   double *x;
   double *y;

   double h=0.00001;
   int n=100000; //How many steps?
   double maxtime = h*double(n);


   vx= new double[n];
   vy= new double[n];
   x= new double[n];
   y= new double[n];

   vx[0]=0;
   vy[0]=2*pi;
   x[0]=1;
   y[0]=0;

   //RK4(vx,vy,x,y,n,h);

   //write (x,y,n,"xout.dat");

   planet erde(10e-6,17,3,4,5,n);

    return 0;
}

void RK4(double *vx,double *vy,double *x,double *y,int n,double h)
{
    double k1vx,k2vx,k3vx,k4vx,k1vy,k2vy,k3vy,k4vy;
    double k1x,k2x,k3x,k4x,k1y,k2y,k3y,k4y;
    for(int i=0;i<n-1;i++)
    {
        k1vx=h*forcex(x[i],y[i]);
        k1vy=h*forcey(x[i],y[i]);
        k2vx=h*forcex(x[i]+k1vx/2.,y[i]+k1vy/2.);
        k2vy=h*forcey(x[i]+k1vx/2.,y[i]+k1vy/2.);
        k3vx=h*forcex(x[i]+k2vx/2.,y[i]+k2vy/2.);
        k3vy=h*forcey(x[i]+k2vx/2.,y[i]+k2vy/2.);
        k4vx=h*forcex(x[i]+k3vx,y[i]+k3vy);
        k4vy=h*forcey(x[i]+k3vx,y[i]+k3vy);

       vx[i+1]=vx[i]+1./6.*(k1vx+2*k2vx+2*k3vx+k4vx);
       vy[i+1]=vy[i]+1./6.*(k1vy+2*k2vy+2*k3vy+k4vy);

        k1x=h*vx[i];
        k2x=h*(vx[i]+k1x/2.);
        k3x=h*(vx[i]+k2x/2.);
        k4x=h*(vx[i]+k3x);
        k1y=h*vy[i];
        k2y=h*(vy[i]+k1y/2.);
        k3y=h*(vy[i]+k2y/2.);
        k4y=h*(vy[i]+k3y);

       x[i+1]=x[i]+1./6.*(k1x+2*k2x+2*k3x+k4x);
       y[i+1]=y[i]+1./6.*(k1y+2*k2y+2*k3y+k4y);

    }
return;
}

double forcex(double x,double y)
{
    double m_S=1;
    double r=sqrt(x*x+y*y);
    return -G*m_S*x/(r*r*r);
}

double forcey(double x,double y)
{
    double m_S=1;
    double r=sqrt(x*x+y*y);
    return -G*m_S*y/(r*r*r);
}
//writes two arrays in a .dat file

void write(double *z, double *y, int n, char *file)
{
    ofstream resout;
    resout.open(file);
    for (int i=0; i<n; i++)
    {
        resout << setprecision(15) << setw(19) << z[i] << " " << setprecision(15) << setw(19) << y[i] << endl;
    }
    resout.close();
}
