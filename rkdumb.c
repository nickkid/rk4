#include "nrutil.h"

//double **y, *xx;//For communication back to main.

/*
Starting from initial values vstart[1..nvar] known 
at x1 use fourth-order Runge-Kutta to advance nstep 
equal increments to x2. The user-supplied routine 
derivs(x,v,dvdx) evaluates derivatives. Results are 
stored in the global variables y[1..nvar][1..nstep+1] 
and xx[1..nstep+1].
*/

void rkdumb(double vstart[], int nvar, double x1, double x2, int nstep,
    double **y, double *xx, void (*derivs)(double, double [], double []))
{
    #include "rk4.h"
    int i, k;
    double x, h;
    double *v, *vout, *dv;

    #ifdef DEBUG
    
    #endif

    v = vector(1, nvar);
    vout = vector(1, nvar);
    dv = vector(1, nvar);
    for (i=1; i<=nvar; i++)
    {
        v[i] = vstart[i];
        y[i][1] = v[i];
    }

    xx[1] = x1;
    x = x1;
    h = (x2 - x1)/nstep;
    for (k=1; k<=nstep; k++)
    {
        (*derivs)(x, v, dv);
        rk4(v, dv, nvar, x, h, vout, derivs);
        if ((double)(x + h) == x) nrerror("Step size too small in routine rkdumb");
        x += h;
        xx[k+1] = x;
        for (i=1; i<=nvar; i++)
        {
            v[i] = vout[i];
            y[i][k+1] = v[i];
        }
    }

    free_vector(dv, 1 ,nvar);
    free_vector(vout, 1, nvar);
    free_vector(v, 1, nvar);
}