extern double k1, k2, k3, k4, k5;
void derivs(double t, double f[], double dfdt[])
{
    /*
    f[1]: concentration of NO
    f[2]: concentration of sGC
    f[3]: concentration of 6C-sGC-NO
    f[4]: concentration of 5C-sGC-NO

    dfdt[1] = -(k1*[NO][sGC] - k2*[6C-sGC-NO])
    dfdt[2] = -(k1*[NO][sGC] - k2*[6C-sGC-NO])
    dfdt[3] = -(dfdt[1] + k3*[NO]*[6C-sGC-NO] - k4*[5C-sGC-NO] + k5*[6C-sGC-NO])
    dfdt[4] = -(k4*[5C-sGC-NO]- k3*[NO][6C-sGC-NO] - k5*[6C-sGC-NO])
    */
   /*
   dfdt[1] = dfdt[2] = k2*f[3]-k1*f[1]*f[2];
   dfdt[3] = -dfdt[1] - k3*f[1]*f[3] + k4*f[4] - k5*f[3];
   dfdt[4] = k3*f[1]*f[3] + k5*f[3] - k4*f[4];
   */

   #ifdef DEBUG
   printf("\nIn dervis.\n");
   #endif
   double r1, r2, r3, r4, r5;
   r1 = k1*f[1]*f[2];
   r2 = k2*f[3];
   r3 = k3*f[1]*f[3];
   r4 = k4*f[4];
   r5 = k5*f[3];
   #ifdef DEBUG
   printf("r1-r5 completed.\n");
   #endif
   dfdt[1] = -r1 + r2;
   dfdt[2] = -r1 + r2;
   dfdt[3] = r1 - r2 - r3 + r4 - r5;
   dfdt[4] = r3 - r4 + r5;
   #ifdef DEBUG
   printf("%f %f %f %f\n", dfdt[1], dfdt[2], dfdt[3], dfdt[4]);
   #endif
}