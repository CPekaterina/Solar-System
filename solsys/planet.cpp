#include "planet.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

planet::planet()
{
    mass =1.0;
    R=new vector[1];
    V=new vector[1];
    R[0].x=0;
    R[0].y=0;
    V[0].x=0;
    V[0].y=0;
}

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
planet::planet(double m, double r,double v, double theta, int n)
{
    mass =m;
    R=new vector[n];
    V=new vector[n];
    R[0].x=r*cos(theta);
    R[0].y=r*sin(theta);
    V[0].x=v*sin(theta);
    V[0].y=-v*cos(theta);
}
void planet::com(planet *planets,int p)
{
    vector P{0,0};//momentum
    vector C{0,0};//mass
    for (int i=1;i<p;i++)//sums up all the moments
    {
        P.x+=planets[i].mass*planets[i].V[0].x;
        P.y+=planets[i].mass*planets[i].V[0].y;

        C.x+=planets[i].mass*planets[i].R[0].x;
        C.y+=planets[i].mass*planets[i].R[0].y;
    }
    R[0].x=-C.x/mass;
    R[0].y=-C.y/mass;
    V[0].x=-P.x/mass; //set the momentum for the first object
    V[0].y=-P.y/mass;

}

void planet::RK4(planet *planets, int p, int n, double h)
{

    for(int i=0;i<n-1;i++) //loop other the number of steps
    {
        for(int j=0; j<p;j++) //loop over the number of objects
        {
            vector f{0.0,0.0}; //intialize force vector
            RK4step(&f,planets,p,j,h,i); //compute the next step for given planet
        }
    }

}

struct vector {double x; double y;};

void planet::force(planet *planets,int p, int j, int i,vector *f,vector l)
{
    f->x=0;//intialize force
    f->y=0;
    double x= planets[j].R[i].x+l.x; //position of considered planet with possible RK4step correction l
    double y= planets[j].R[i].y+l.y;
    for(int k=0; (k<p);k++) //loop over all planets except for the considered one
    {
        if(k!=j) //excludes the considered planet
        {
            double dx=x-planets[k].R[i].x; //distance to planet k
            double dy=y-planets[k].R[i].y;
            double r=sqrt(dx*dx+dy*dy);

            f->x += -G*planets[k].mass*dx/(r*r*r); //gravitational force
            f->y += -G*planets[k].mass*dy/(r*r*r);

        }
    }
}
void planet::RK4step(vector *f,planet *planets, int p, int j,double h, int i)
{
    vector k1v,k2v,k3v,k4v;
    vector k1,k2,k3,k4;
    vector l{0.0,0.0};

    double Vx=planets[j].V[i].x;
    double Vy=planets[j].V[i].y;

    force(planets,p,j,i,f,l);

    k1v.x=h/2.*f->x;
    k1v.y=h/2.*f->y;

    k1.x=h/2.*Vx;
    k1.y=h/2.*Vy;

    force (planets, p,j,i,f,k1);

    k2v.x=h/2.*f->x;
    k2v.y=h/2.*f->y;
    k2.x=h/2.*(Vx+k1v.x);
    k2.y=h/2.*(Vy+k1v.y);

    force(planets, p,j,i,f,k2);

    k3v.x=h*f->x;
    k3v.y=h*f->y;
    k3.x=h*(Vx+k2v.x);
    k3.y=h*(Vy+k2v.y);

    force(planets, p,j,i,f,k3);

    k4v.x=h*f->x;
    k4v.y=h*f->y;
    k4.x=h*(Vx+k3v.x);
    k4.y=h*(Vy+k3v.y);

   //write the next positions and velocities

   planets[j].V[i+1].x=planets[j].V[i].x+1./6.*(2.*k1v.x+4*k2v.x+2.*k3v.x+k4v.x);
   planets[j].V[i+1].y=planets[j].V[i].y+1./6.*(2.*k1v.y+4*k2v.y+2.*k3v.y+k4v.y);


   planets[j].R[i+1].x=planets[j].R[i].x+1./6.*(2.*k1.x+4.*k2.x+2.*k3.x+k4.x);
   planets[j].R[i+1].y=planets[j].R[i].y+1./6.*(2.*k1.y+4.*k2.y+2.*k3.y+k4.y);

}

void planet::RXYwrite(int n, char *file)
{
    ofstream resout;
    resout.open(file);
    for (int i=0; i<n; i++)
    {
        resout << setprecision(15) << setw(19) << R[i].x << " " << setprecision(15) << setw(19) << R[i].y << endl;
    }
    resout.close();
}


void planet::VXYwrite(int n, char *file)
{
    ofstream resout;
    resout.open(file);
    for (int i=0; i<n; i++)
    {
        resout << setprecision(15) << setw(19) << V[i].x << " " << setprecision(15) << setw(19) << V[i].y << endl;
    }
    resout.close();

}


void planet::Verlet(planet *planets, int p, int n, double h)
{
    vector l{0,0};

    for(int j=0;j<p;j++) //not self starting -> the second step is computed with taylor expansion
    {
        vector f;

        force(planets,p,j,0,&f,l);

        planets[j].R[1].x = planets[j].R[0].x + planets[j].V[0].x*h + 0.5*h*h*f.x;
        planets[j].R[1].y = planets[j].R[0].y + planets[j].V[0].y*h + 0.5*h*h*f.y;
    }

    for(int i=1;i<n-1;i++) //Verlet iteration for all time steps
    {
        for(int j=0;j<p;j++){ //and all planets

            vector f;
            force(planets,p,j,i,&f,l);

            planets[j].R[i+1].x = 2.*planets[j].R[i].x - planets[j].R[i-1].x + h*h*f.x;
            planets[j].R[i+1].y = 2.*planets[j].R[i].y - planets[j].R[i-1].y + h*h*f.y;

            planets[j].V[i].x=(planets[j].R[i+1].x-planets[j].R[i-1].x)/(2*h); //velocities from the 2point formula
            planets[j].V[i].y=(planets[j].R[i+1].y-planets[j].R[i-1].y)/(2*h);
        }
    }

}



