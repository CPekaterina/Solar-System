#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <planet.h>

using namespace std;
const double pi= 3.14159265359;
const double G=4.*pi*pi;
<<<<<<< HEAD
=======

double forcex(double x,double y);
double forcey(double x,double y);
>>>>>>> origin/master
void write(double *z, double *y, int n, char *file);
//time - year, mass - sun mass, distance - AU
<<<<<<< HEAD

int main()
{

   double h=1e-3;
   int n=1200e3; //How many steps?
   double maxtime = h*double(n);

   //planet definitions

   planet planets[10];

   planet erde(3*1e-6,1,2*pi,8.7,n);
   planet sonne(1,0,0,0,0,n);
   planet jupiter(1.9/2.*1e-3,5.29,10./12.*pi,127.5,n);
   planet saturn(5.5/2.*1e-4,9.94,2.037,234.5,n);
   planet mars(3.3*1e-7,1.52,5.073,291,n);
   planet venus(4.9/2.*1e-6,0.72,7.36,173,n);
   planet uranus(4.4*1e-5,19.19,1.4317,14,n);
   planet merkur(1.2*1e-7,0.41,9.957,317.9,n);
   planet neptun(1.03/2.*1e-7,30.06,1.142,336,n);
   planet pluto(1.31/2.*1e-8,32.73,0.9923,282.7,n);
   planets[0]=sonne;
   planets[1]=merkur;
   planets[2]=venus;
   planets[3]=erde;
   planets[4]=mars;
   planets[5]=jupiter;
   planets[6]=saturn;
   planets[7]=uranus;
   planets[8]=neptun;
   planets[9]=pluto;



sonne.com(planets,10);
cout << sonne.V[0].x << " " << sonne.V[0].y << endl;




 erde.RK4(planets,10,n,h);
// erde.Verlet(planets,2,n,h);


// Energy conservation
/*

   double *en,*time;
   en= new double [n];
   time = new double[n];
   for (int i=0;i<n;i++)
   {
       double qr=erde.R[i].x*erde.R[i].x+erde.R[i].y*erde.R[i].y;
       en[i]=0.5*3*1e-6*(erde.V[i].x*erde.V[i].x+erde.V[i].y*erde.V[i].y); //kin
       //en[i]=G*3*1e-6/sqrt(qr);                                             //pot
       //en[i]=3*1e-6*(erde.R[i].x*erde.V[i].y-erde.R[i].y*erde.V[i].x);
       time[i]=i*h;
   }
   write(time,en,n,"kinenergV.dat");
*/

//data production


   erde.RXYwrite(n,"3.dat");
   sonne.RXYwrite(n,"0.dat");
    jupiter.RXYwrite(n,"5.dat");
    merkur.RXYwrite(n,"1.dat");
    venus.RXYwrite(n,"2.dat");
    saturn.RXYwrite(n,"6.dat");
    uranus.RXYwrite(n,"7.dat");
    neptun.RXYwrite(n,"8.dat");
    pluto.RXYwrite(n,"9.dat");
    mars.RXYwrite(n,"4.dat");


=======
void RK4(double *vx,double *vy,double *x,double *y,int n,double h);
void Verlet(double *vx, double *vy, double *x, double *y, int n, double h);

int main()
{
    double m_S=1.;
    double *vx;
    double *vy;
    double *x;
    double *y;
    double h=1e-4;
    int n=5e4; //How many steps?
    double maxtime = h*double(n);

    //planet definitions
    planet planets[3];
    planet erde(3*1e-6,1.,0,0,2.*pi,n);
    planet sonne(1,0,0,0,0,n);
    planet jupiter(1,5,0,0,10./12.*pi,n);
    planets[0]=erde;
    planets[1]=sonne;
    planets[2]=jupiter;
    // erde.RK4(planets,2,n,h);
    // erde.Verlet(planets,2,n,h);
    // Energy conservation
    double *en,*time;
    en= new double [n];
    time = new double[n];
    for (int i=0;i<n;i++)
    {
        double qr=erde.R[i].x*erde.R[i].x+erde.R[i].y*erde.R[i].y;
        //en[i]=0.5*3e-6*(erde.V[i].x*erde.V[i].x+erde.V[i].y*erde.V[i].y); //kin
        //en[i]=1/qr; //pot
        time[i]=i*h;
    }
    write(time,en,n,"kinenergy.dat");
    //data production
    erde.RXYwrite(n,"erdeV.dat");
    // sonne.RXYwrite(n,"sonne.dat");
    // jupiter.RXYwrite(n,"jupiter.dat");
>>>>>>> origin/master
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
