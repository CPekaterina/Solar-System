#ifndef PLANET_H
#define PLANET_H


class planet
{
private:
    double mass;
    double pi = 3.14159265359;
    double G = 4.*pi*pi;
    double c= 63072000;

public:
    planet();                                                           //default constructor
    planet(double m, int n);                                            //c. including mass and number of steps
    planet(double m, double r,double v, double theta, int n);           //c. including positions and velocity in polar coordinates
    planet(double m, double x, double y, int n);                        //c. including cartesian positions
    planet(double m, double x, double y, double vx, double vy, int n);  //c. including cartesian positions and velocities

    struct vector {double x; double y;};

    vector *R;
    vector *V;
    vector f;

    void RK4(planet *planets, int p, int n, double h);
    /* Receives an array of the classes' objects (*planets) and performs RK4step for a given number of time steps n
     *  of certain step length h and the number of objects p which is also the size of the array.
     * Then the velocities and positions for the next time step are computed for all planets before moving on to the next time step.
    */

    void RK4step(vector *f,planet *planets, int p, int j,double h, int i);
    /* Uses the values for the current time step i and planet j to compute the velocity
     * and position for next time step for said planet taking in account the forces that the other p-1 objects stored in *planets
     * generate . Here the RK4 algorithm is applied and the function force is used.
     */


    void force(planet *planets,int p, int j, int i,vector *f,vector l);
    /* Computes the gravitational force and stores it on *f for a certain planet j experiences given the current positions
     * and masses of the other p objects in *planets involved. l is a possible correction to the position due to the RK4 slope
     * computation
     */

    void RXYwrite(int n, char *file);
    void VXYwrite(int n, char *file);
    /* write the positions or the velocities for n timepsteps into a .dat file named *file.
     */

    void Verlet(planet *planets, int p, int n, double h);
    /* Receives the same data as RK4 but performs the method directly inside the function and also computes the velocities.
     */

    void write(double *z, double *y, int n, char *file);
    /*general write function
     */

    void com(planet *planets,int p);
    /* Computes the initial velocity for the first object so that the total momentum of the system is zero.
     */

};

#endif // PLANET_H
