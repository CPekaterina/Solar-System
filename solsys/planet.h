#ifndef PLANET_H
#define PLANET_H


class planet
{
private:
    double mass;
    double pi = 3.1425926;
    double G = 4*pi*pi;

public:
    planet(double m, int n);
    planet(double m, double x, double y, int n);
    planet(double m, double x, double y, double vx, double vy, int n);

    struct vector {double x; double y;};

    vector *R;
    vector *V;
    //vector f;

    void RK4(planet *planets, int p, int n, double h);
    void force(planet *planets, int p, int j, int i, vector f);
    void RK4step(vector force,planet *planets, int j, int i, double h);

};

#endif // PLANET_H
