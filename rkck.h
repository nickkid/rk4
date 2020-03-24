#ifndef RKCK_H_INCLUDED
#define RKCK_H_INCLUDED
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
    double yerr[], void (*derivs)(double, double [], double []));
#endif