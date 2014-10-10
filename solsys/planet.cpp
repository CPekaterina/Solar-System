#include "planet.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

planet::planet()
{
    mass =1;
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

void planet::RK4(planet *planets, int p, int n, double h)
{

    for(int i=0;i<n-1;i++)
    {
        for(int j=0; j<p;j++)
        {
            vector f{0,0};
            RK4step(&f,planets,p,j,h,i);
           //cout << "Kraft: "<<f.x << " " << f.y << endl;
        }
    }
}

struct vector {double x; double y;};

void planet::force(planet *planets,int p, int j, int i,vector *f,vector l)
{

//  cout << "p: " << p << ", j: " << j << ", i: " << i<< endl;
//  cout << "l.x: " <<l.x << ", l.y: " << l.y << endl;
    f->x=0;
    f->y=0;
    double x= planets[j].R[i].x+l.x;
    double y= planets[j].R[i].y+l.y;
    for(int k=0; (k<p);k++)
    {
        if(k!=j)
        {
            double dx=x-planets[k].R[i].x;
            double dy=y-planets[k].R[i].y;
          //cout << "dx: " << dx << ", dy:" << dy <<endl;

            double r=sqrt(dx*dx+dy*dy);
            f->x += -G*planets[k].mass*dx/(r*r*r);
            f->y += -G*planets[k].mass*dy/(r*r*r);
          // cout << "zwischenKraft: "<<f->x << " " << f->y << endl;
            //cout << planets[k].mass << endl;
            //cout << planets[j].mass << endl;

        }

     }

}
void planet::RK4step(vector *f,planet *planets, int p, int j,double h, int i)
{
    vector k1v,k2v,k3v,k4v;
    vector k1,k2,k3,k4;

    vector l{0,0};

    force(planets,p,j,i,f,l);

    k1v.x=h/2.*f->x;
    k1v.y=h/2.*f->y;

    //cout << "force_x: " <<f->x << ", force_y: " <<f->y << endl;


    force (planets, p,j,i,f,k1v);

    k2v.x=h/2.*f->x;
    k2v.y=h/2.*f->y;

    //cout << "force_x: " <<f->x << ", force_y: " <<f->y << endl;

    force(planets, p,j,i,f,k2v);

    k3v.x=h*f->x;
    k3v.y=h*f->y;

    //cout << "force_x: " <<f->x << ", force_y: " <<f->y << endl;

    force(planets, p,j,i,f,k3v);

    k4v.x=h*f->x;
    k4v.y=h*f->y;

    //cout << "force_x: " <<f->x << ", force_y: " <<f->y << endl;

   planets[j].V[i+1].x=planets[j].V[i].x+1./6.*(2*k1v.x+4*k2v.x+2*k3v.x+k4v.x);
   planets[j].V[i+1].y=planets[j].V[i].y+1./6.*(2*k1v.y+4*k2v.y+2*k3v.y+k4v.y);

    double Vx=planets[j].V[i].x;
    double Vy=planets[j].V[i].y;
    //cout << "j: " << j <<  "  Vx: " <<Vx << ", Vy: " << Vy <<endl;
    k1.x=h*Vx;
    k2.x=h*(Vx+k1.x/2.);
    k3.x=h*(Vx+k2.x/2.);
    k4.x=h*(Vx+k3.x);
    k1.y=h*Vy;
    k2.y=h*(Vy+k1.y/2.);
    k3.y=h*(Vy+k2.y/2.);
    k4.y=h*(Vy+k3.y);

   planets[j].R[i+1].x=planets[j].R[i].x+1./6.*(k1.x+2*k2.x+2*k3.x+k4.x);
   planets[j].R[i+1].y=planets[j].R[i].y+1./6.*(k1.y+2*k2.y+2*k3.y+k4.y);

   //cout << "j: " << j <<  "  Rx: " <<   planets[j].R[i+1].x << ", Ry: " <<  planets[j].R[i+1].y <<endl;

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

