#ifndef RKDUMB_H_INCLUDE
#define RKDUMB_H_INCLUDE

void rkdumb(double vstart[], int nvar, double x1, double x2, int nstep,
    double **y, double *xx, void (*derivs)(double, double [], double []));
    
#endif