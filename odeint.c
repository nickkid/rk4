#include <math.h>
#include "odeint.h"
#include "nrutil.h"

#ifndef DEBUG
//#define DEBUG
#endif

#ifdef DEBUG
#include <stdio.h>
#endif

#define MAXSTP 10000000
#define TINY 1.0e-30

extern int kmax, kount;
extern double *xp, **yp, dxsav;
/*
User storage for intermediate results. Preset kmax and dxsav in the calling program. If kmax ?= 
0 results are stored at approximate intervals dxsav in the arrays xp[1..kount], yp[1..nvar][1..kount], 
where kount is output by odeint. Defining declarations for these variables, with memory allocations 
xp[1..kmax] and yp[1..nvar][1..kmax] for the arrays, should be in the calling program.
*/

/*
Runge-Kutta driver with adaptive stepsize control. Integrate starting values ystart[1..nvar] 
from x1 to x2 with accuracy eps, storing intermediate results in global variables. h1 should 
be set as a guessed first stepsize, hmin as the minimum allowed stepsize (can be zero). On 
output nok and nbad are the number of good and bad (but retried and fixed) steps taken, and 
ystart is replaced by values at the end of the integration interval. derivs is the user-supplied 
routine for calculating the right-hand side derivative, while rkqs is the name of the stepper 
routine to be used.
*/

void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
    double hmin, int *nok, int *nbad, 
    void (*derivs)(double, double [], double []),
    void (*rkqs)(double [], double [], int, double *, double, double, double [],
    double *, double *, void (*)(double, double [], double [])))

{
    #ifdef DEBUG
    printf("\nIn odeint.\n");
    printf("Variables in odeint:\n");
    printf("nvar = %d\nx1 = %f\nx2 = %f\neps = %g\nh1 = %g\nhmin = %g\n",
    nvar, x1, x2, eps, h1, hmin);
    #endif
    int nstp,i;
    double xsav,x,hnext,hdid,h;
    double *yscal,*y,*dydx;
    
    yscal = vector(1, nvar);
    y = vector(1, nvar);
    dydx = vector(1, nvar);
    x = x1;
    h=SIGN(h1, x2-x1);
    *nok = (*nbad) = kount = 0;
    for (i=1;i<=nvar;i++) y[i] = ystart[i];
    if (kmax > 0) xsav = x - dxsav*2.0;//Assures storage of first step.
    #ifdef DEBUG
        printf("h = %g\n", h);
        printf("Begin iterations in odeint.\n");
        printf("nvar = %d in odeint.\n", nvar);
    #endif
    //for (nstp=1; nstp <= MAXSTP; nstp++)
    for (nstp=1;; nstp++)
    {
        #ifdef DEBUG
            printf("nstp = %d.\n", nstp);
            for (int i=1; i<=nvar; i++)
                printf("y[%d] = %g\n", i, y[i]);
            #endif
        for (i=1; i <= nvar; i++)
        //Scaling used to monitor accuracy. This general-purpose choice can be modified if need be.
            yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + TINY;
        if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav))
        {
            xp[++kount] = x;//Store intermediate results.
            for (i=1; i <= nvar; i++) yp[i][kount]=y[i];
            xsav=x;
        }
        if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;//If stepsize can overshoot, decrease.
        (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
        #ifdef DEBUG
            printf("rkqs in odeint iteration %d completed.\n", nstp);
        #endif
        if (hdid == h) ++(*nok); else ++(*nbad);
        if ((x-x2)*(x2-x1) >= 0.0)//Are we done?
        {
            for (i=1; i <= nvar; i++) ystart[i]=y[i];
            if (kmax)
            {
                xp[++kount] = x;//Save final step.
                for (i=1; i <= nvar; i++) yp[i][kount]=y[i];
            }
            free_vector(dydx,1,nvar);
            free_vector(y,1,nvar);
            free_vector(yscal,1,nvar);
            return;//Normal exit.
        }
        #ifdef DEBUG
        printf("hnext = %g\nhmin = %g\n", hnext, hmin);
        #endif
        if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
        h=hnext;
    }
    nrerror("Too many steps in routine odeint");
}