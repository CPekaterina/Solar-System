#ifndef PLANET_H
#define PLANET_H


class planet
{
private:
    double mass;
    double pi = 3.14159265359;
    double G = 4.*pi*pi;

public:
    planet();
    planet(double m, int n);
    planet(double m, double r,double v, double theta, int n);
    planet(double m, double x, double y, int n);
    planet(double m, double x, double y, double vx, double vy, int n);

    struct vector {double x; double y;};

    vector *R;
    vector *V;
    vector f;

    void RK4(planet *planets, int p, int n, double h);
    void force(planet *planets,int p, int j, int i,vector *f,vector l);
    void RK4step(vector *f,planet *planets, int p, int j,double h, int i);
    void RXYwrite(int n, char *file);
    void VXYwrite(int n, char *file);
    void Verlet(planet *planets, int p, int n, double h);
    void write(double *z, double *y, int n, char *file);
    void com(planet *planets,int p);

};

#endif // PLANET_H
