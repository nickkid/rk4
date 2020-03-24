#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "rkqs.h"
#include "odeint.h"
#include "derivs.h"
//#define DEBUG

//#define DEBUG
/*dimensions specification:
time: s
concentration for k's: M
concentration for sGC: uM
concentration for NO: uM
k1, k3: M^-1
k2, k4, k5: s^-1 
epsilon: mM^-1*cm^-1
*/

const int nvar = 4;//NO, sGC, 6C, 5C
int nok, nbad, kount;
int kmax = 10000000;
int current_case = 0;

double **yp, *xp;
double *vstart;
double x1, x2;
double dxsav =  1e-4;
/*
const double k1 = 1.55e8;
const double k2 = 1e-2;
const double k3 = 3.3e5;
const double k4 = 1e-3;
const double k5 = 5e-4;
const double length = 1;
*/
double k1, k2, k3, k4, k5, length;
const double epsilons[] = {124e3, 66e3, 46e3};
void init(char* input_file)
{
    nok = 0;
    nbad = 0;
    kount = 0;
    kmax = 10000000;
    yp = matrix(1, nvar, 1, kmax);

    #ifdef DEBUG
    printf("yp created.\n");
    #endif

    vstart = vector(1, nvar);

    #ifdef DEBUG
    printf("vstart created.\n");
    #endif
    FILE *file = fopen(input_file, "r");
    printf("%s\n", input_file);
    //initial concentrations
    fscanf(file, "%lf", &k1);
    fscanf(file, "%lf", &k2);
    fscanf(file, "%lf", &k3);
    fscanf(file, "%lf", &k4);
    fscanf(file, "%lf", &k5);
    fscanf(file, "%lf", &length);
    fscanf(file, "%lf", &vstart[1]); //vstart[1] = [NO]
    fscanf(file, "%lf", &vstart[2]); //vstart[2] = [sGC]
    fscanf(file, "%lf", &vstart[3]); //vstart[3] = [6C-sGC-NO]
    fscanf(file, "%lf", &vstart[4]); //vstart[4] = [5C-sGC-NO]
    //start and end time of simulation
    fscanf(file, "%lf", &x1);
    fscanf(file, "%lf", &x2);

    printf("k1 = %lf\n", k1);
    printf("k2 = %lf\n", k2);
    printf("k3 = %lf\n", k3);
    printf("k4 = %lf\n", k4);
    printf("k5 = %lf\n", k5);
    printf("length = %lf\n", length);
    printf("[NO] = %lf\n", vstart[1]);
    printf("[sGC] = %lf\n", vstart[2]);
    printf("[6C-sGC-NO] = %lf\n", vstart[3]);
    printf("[5C-sGC-NO] = %lf\n", vstart[4]);
    printf("t0 = %lf\n", x1);
    printf("t1 = %lf\n", x2);
    xp = vector(1, kmax);
    #ifdef DEBUG
    printf("xp created.\n");
    #endif
    fclose(file);
}

void save_results(char* postfix)
{
    char const *prefix = "results";
    char const *extension = ".csv";
    static const int filename_maximum_length = 100;
    /*
    char string_icase[8] = "0";
    itoa(icase, string_icase, 10);
    */
    char filename[filename_maximum_length];
    sprintf(filename, "%s%s%s\0", prefix, postfix, extension);
    FILE *file = fopen(filename, "w");
    
    if (file == NULL) return;
    #ifdef DEBUG
    printf("Results file created!\n");
    #endif
    /*
    for (int i=0; i < ncase; i++)
        for (int j=1; j <= nstep + 1; j++)
        fprintf(file, "%f,%f,%f,%f,%f\n", xx[j], y[i][1][j], y[i][2][j], y[i][3][j], y[i][4][j]);
    */
    fprintf(file, "Time,[NO],[sGC],[6C-sGC-NO],[5C-sGC-NO],Absorbance\n");
    for (int istep=1; istep <= kount; istep ++)
    {
        double Absorbance = 0;
        fprintf(file, "%e", xp[istep]);
        for (int ivar=1; ivar <= nvar; ivar ++)
        {
            fprintf(file, ",%e", yp[ivar][istep]);
            if (ivar != 1)
                Absorbance += epsilons[ivar-2] * yp[ivar][istep];
        }
        Absorbance *= length;
        fprintf(file, ",%e\n", Absorbance);   
    }
    fclose(file);
    return ;
}

void free_space()
{
    free_matrix(yp, 1, nvar, 1, kmax);
    free_vector(xp, 1, kmax);
}
/*
argv[]:
0: name of executable
1: input file name
2: output file order
*/
int main(int argc, char *argv[])
{
    printf("%d ", argc);
    for (int i=0; i<argc; i++)
        printf("%s ", argv[i]);
    #ifdef DEBUG
    printf("In main.\n");
    #endif
    printf("strlen() = %d ", strlen(argv[1]));
    init(argv[1]);
    double eps = 1e-3;
    double h1 = 1e-3;
    double hmin = 1e-40;
    #ifdef DEBUG
        printf("k1 = %lf\nk2 = %lf\nk3 = %lf\nk4 = %lf\nk5 = %lf\nlength = %lf\n",
        k1, k2, k3, k4, k5, length);
        printf("nvar = %d\nx1 = %f\nx2 = %f\neps = %g\nh1 = %g\nhmin = %g\n",
    nvar, x1, x2, eps, h1, hmin);
        printf("\n");
    #endif
    odeint(vstart, nvar, x1, x2, eps, h1, hmin, &nok, &nbad,
        derivs,
        rkqs);
    #ifdef DEBUG
        printf("odeint completed.\n");
    #endif
    save_results(argv[2]);
    free_space();
    
    #ifdef DEBUG
        printf("Calculation complete!\n");
    #endif
    return 0;
}
