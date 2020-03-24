#include <math.h>
#include "rkqs.h"
#include "rkck.h"
#include "nrutil.h"
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.90e-4

#ifndef DEBUG
//#define DEBUG
#endif

#ifdef DEBUG
#include <stdio.h>
#endif

/*
Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and adjust stepsize. 
Input are the dependent variable vector y[1..n] and its derivative dydx[1..n] at the starting value of the 
independent variable x. Also input are the stepsize to be attempted htry, the required accuracy eps, and 
the vector yscal[1..n] against which the error is scaled. On output, y and x are replaced bytheir new values,
 hdid is the stepsize that was actually accomplished, and hnext is the estimated next stepsize. derivs is the 
 user-supplied routine that computes the right-hand side derivatives.
*/
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
    double yscal[], double *hdid, double *hnext,
    void (*derivs)(double, double [], double []))
{
    /*
    The following declaration gets me stuck for a long time. Variable type of n is wrong.
    void rkck(double y[], double dydx[], double n, double x, double h,
        double yout[], double yerr[], void (*derivs)(double, double [], double []));
    */
   /*
    void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
        double yerr[], void (*derivs)(double, double [], double []));
        */
    int i;
    double errmax, h, htemp, xnew, *yerr, *ytemp;

    yerr = vector(1, n);
    ytemp = vector(1, n);
    h = htry;                                           //Set stepsize to the initial trial value
    #ifdef DEBUG
    printf("\nIn rkqs.\n");
    printf("htry = %g\n", htry);
    printf("eps = %g\n", eps);
    int it = 0;
    #endif
    for (;;)
    {
        #ifdef DEBUG
            it ++;
            printf("\nIteration %d in rkqs begins.\n", it);
            printf("nvar = %d.\n", n);
            printf("x = %g\n", *x);
        #endif

        rkck(y, dydx, n, *x, h, ytemp, yerr, derivs);   //Take a step.
        #ifdef DEBUG
            printf("h = %g\n", h);
        #endif
        #ifdef DEBUG
            printf("rkck in rkqs completed.\n");
        #endif
        
        errmax = 0.0;                                   //Evaluate accuracy
        for (i=1; i<=n; i++)
        {
            errmax = DMAX(errmax, fabs(yerr[i]/yscal[i]));
            #ifdef DEBUG
                printf("yerr[%d] = %g   yscal[%d] = %g\t\
                 errmax = %g\n", i, yerr[i], i, yscal[i], errmax);
            #endif
        }
        
        #ifdef DEBUG
            printf("errmax = %g\n", errmax);
        #endif
        
        errmax /= eps;                                  //Scale relative to required tolerance.
        #ifdef DEBUG
            printf("errmax = %g\n", errmax);
        #endif
        if (errmax <= 1.0) break;                       //Step succeeded. Compute size of next step.
        htemp = SAFETY*h*pow(errmax, PSHRNK);
        h = (h >= 0.0 ? DMAX(htemp, 0.1*h) : DMIN(htemp, 0.1*h));
                                                        //No more than a factor of 10.
        xnew = (*x) + h;
        if (xnew == *x) nrerror("stepsize underflow in rkqs");
    }

    if (errmax > ERRCON) *hnext = SAFETY*h*pow(errmax, PGROW);
    else *hnext = 5.0*h;                                //No more than a factor 5 increase.
    *x += (*hdid=h);
    for (i=1; i<=n; i++) y[i]=ytemp[i];
    free_vector(ytemp, 1, n);
    free_vector(yerr, 1, n);
}