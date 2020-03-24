#include "nrutil.h"

/*
Given values for the variables y[1..n] and their derivatives dydx[1..n] 
known at x,use the fourth-order Runge-Kutta method to advance the 
solution over an interval h and return the incremented variables as 
yout[1..n], which need not be a distinct array from y. The user supplies 
the routine derivs(x,y,dydx), which returns derivatives dydx at x.
*/
void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
    void (*derivs)(double, double [], double[]))
{
    /*
    k1 = h*f(xn, yn)
    k2 = h*f(xn + h/2, yn + k1/2)
    k3 = h*f(xn + h/2, yn + k2/2)
    k4 = h*f(xn  + h, yn + k3)
    y(n+1) = y(n) + k1/6 + k2/3 + k3/3 + k4/6 + O(h^5)
    */
    int i;
    double xh, hh, h6, *dym, *dyt, *yt;

    dym = vector(1, n);
    dyt = vector(1, n);
    yt = vector(1, n);
    hh = h*0.5;
    h6 = h/6.0;
    xh = x + hh;

    for (i=1; i<=n; i++) yt[i] = y[i] + hh*dydx[i];//First step.

    (*derivs)(xh, yt, dyt);                        //Second step.
    for (i=1; i<=n; i++) yt[i] = y[i] + hh*dyt[i]; 

    (*derivs)(xh, yt, dym);                        //Third step.
    for (i=1; i<=n; i++)
    {
        yt[i] = y[i] + h*dym[i];
        dym[i] += dyt[i];
    }
    
    (*derivs)(x + h, yt, dyt);                    //Fourth step.

    for (i=1; i<=n; i++)                         //Accumulate increments
        yout[i] = y[i] + h6*(dydx[i] + dyt[i] + 2.0*dym[i]);//with proper weights

    free_vector(yt, 1, n);
    free_vector(dyt, 1, n);
    free_vector(dym, 1, n);
}
