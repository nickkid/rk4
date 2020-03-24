#ifndef RKQS_H_INCLUDED
#define RKQS_H_INCLUDED
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
    double yscal[], double *hdid, double *hnext,
    void (*derivs)(double, double [], double []));
#endif