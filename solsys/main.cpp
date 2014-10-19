#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <planet.h>
using namespace std;
const double pi= 3.14159265359;
const double G=4.*pi*pi;
double forcex(double x,double y);
double forcey(double x,double y);
void write(double *z, double *y, int n, char *file);
//time - year, mass - sun mass, distance - AU
void RK4(double *vx,double *vy,double *x,double *y,int n,double h);
void Verlet(double *vx, double *vy, double *x, double *y, int n, double h);
int main()
{
double m_S=1.;
double *vx;
double *vy;
double *x;
double *y;
double h=0.17;
int n=5883; //How many steps?
double maxtime = h*double(n);
//planet definitions
planet planets[3];
planet erde(3*1e-6,1.,0,0,2*pi,n);
planet sonne(1,0,0,0,0,n);
planet jupiter(1e-3,5.2,0,0,2.7478368,n);
planets[0]=erde;
planets[1]=sonne;
planets[2]=jupiter;

//erde.RK4(planets,2,n,h);
//erde.Verlet(planets,2,n,h);

//erde.RK4(planets,3,n,h);
erde.Verlet(planets,3,n,h);
//jupiter.RK4(planets,3,n,h);
jupiter.Verlet(planets,3,n,h);


// Energy conservation
double *en,*time;
en= new double [n];
time = new double[n];
for (int i=0;i<n;i++)
{
double qr=erde.R[i].x*erde.R[i].x+erde.R[i].y*erde.R[i].y;
en[i]=0.5*3e-6*(erde.V[i].x*erde.V[i].x+erde.V[i].y*erde.V[i].y); //kin
//en[i]=1./sqrt(qr); //pot
time[i]=i*h;
}
//distance
double *re;
double *rj;
re = new double[n];
rj = new double[n];

for(int i=0;i<n;i++){
    re[i]=sqrt(erde.R[i].x*erde.R[i].x+erde.R[i].y*erde.R[i].y);
    rj[i]=sqrt(jupiter.R[i].x*jupiter.R[i].x+jupiter.R[i].y*jupiter.R[i].y);
}

write(time,re,n,"distance_earthV.dat");
write(time,rj,n,"distance_jupiterV.dat");

write(time,en,n,"kinenergy.dat");
//data production
erde.RXYwrite(n,"erdeV.dat");
// sonne.RXYwrite(n,"sonne.dat");
jupiter.RXYwrite(n,"jupiterV.dat");
return 0;
}
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
