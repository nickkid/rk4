#ifndef RK4_H_INCLUDED
#define RK4_H_INCLUDED

void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
    void (*derivs)(double, double [], double[]));
    
#endif