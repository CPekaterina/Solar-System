#include "planet.h"
#include <cmath>

planet::planet(double m, int n)
{
    mass =m;
    R=new vector[n];
    V=new vector[n];
    R[0].x=1;
    R[0].y=0;
    V[0].x=0;
    V[0].y=2*pi;
}
planet::planet(double m,double x, double y, int n)
{
    mass =m;
    R=new vector[n];
    V=new vector[n];
    R[0].x=x;
    R[0].y=y;
    V[0].x=0;
    V[0].y=2*pi;
}
planet::planet(double m,double x, double y, double vx, double vy, int n)
{
    mass =m;
    R=new vector[n];
    V=new vector[n];
    R[0].x=x;
    R[0].y=y;
    V[0].x=vx;
    V[0].y=vy;
}

void planet::RK4(planet *planets, int p, int n, double h)
{

    for(int i=0;i<n-1;i++)
    {
        for(int j=0; j<p;j++)
        {
            vector f;
            force(planets,p,j,i,f);
            //RK4step(f,planets,j,h,i);
        }
    }
}

struct vector {double x; double y;};

void planet::force(planet *planets,int p, int j, int i,vector f)
{

    double x= planets[j].R[i].x;
    double y= planets[j].R[i].y;
    for(int k=0; k<p && k!=j;k++)
    {
        double dx=x-planets[k].R[i].x;
        double dy=y-planets[k].R[i].y;

        double r=sqrt(dx*dx+dy*dy);
        f.x += -G*planets[k].mass*dx/(r*r*r);
        f.y += -G*planets[k].mass*dy/(r*r*r);
    }

}
