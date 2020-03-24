#ifndef ODEINT_H_INCLUDED
#define ODEINT_H_INCLUDED
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1, double hmin, int *nok, int *nbad, 
    void (*derivs)(double, double [], double []),
    void (*rkqs)(double [], double [], int, double *, double, double, double [],
    double *, double *, void (*)(double, double [], double [])));
#endif