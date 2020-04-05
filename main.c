#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "nrutil.h"
#include "rkqs.h"
#include "odeint.h"
#include "derivs.h"
#include "loadcsv.h"
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

int number_of_samples = 0;
double **data = NULL;

double **yp, *xp;
double *vstart;
double x1, x2;
double dxsav =  1e-4;
double eps = 1e-3;
double h1 = 1e-3;
double hmin = 1e-40;
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
void init(char*);
void init_parameter(char* []);
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
    fscanf(file, "%lf", &eps);
    fscanf(file, "%lf", &h1);
    fscanf(file, "%lf", &hmin);

    xp = vector(1, kmax);
    #ifdef DEBUG
    printf("xp created.\n");
    #endif
    fclose(file);
}

void init_parameter(char* argv[])
{
    nok = 0;
    nbad = 0;
    kount = 0;
    kmax = 10000000;
    yp = matrix(1, nvar, 1, kmax);
    vstart = vector(1, nvar);

    //initial concentrations
    sscanf(argv[1], "%lf", &k1);
    sscanf(argv[2], "%lf", &k2);
    sscanf(argv[3], "%lf", &k3);
    sscanf(argv[4], "%lf", &k4);
    sscanf(argv[5], "%lf", &k5);
    sscanf(argv[6], "%lf", &length);
    sscanf(argv[7], "%lf", &vstart[1]); //vstart[1] = [NO]
    sscanf(argv[8], "%lf", &vstart[2]); //vstart[2] = [sGC]
    sscanf(argv[9], "%lf", &vstart[3]); //vstart[3] = [6C-sGC-NO]
    sscanf(argv[10], "%lf", &vstart[4]); //vstart[4] = [5C-sGC-NO]
    //start and end time of simulation
    sscanf(argv[11], "%lf", &x1);
    sscanf(argv[12], "%lf", &x2);
    sscanf(argv[13], "%lf", &eps);
    sscanf(argv[14], "%lf", &h1);
    sscanf(argv[15], "%lf", &hmin);

    xp = vector(1, kmax);
    #ifdef DEBUG
    printf("xp created.\n");
    #endif
}

void save_results(char* filename)
{
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

double calculate_std()
{
    double std = 0;
    for (int istep=1; istep <= kount; istep ++)
    {
        double Absorbance = 0;
        for (int ivar=2; ivar <= nvar; ivar ++)
            Absorbance += epsilons[ivar-2] * yp[ivar][istep];
        Absorbance *= length;
        double difference = Absorbance - data[istep][2];
        std += (Absorbance * Absorbance);
    }
    std = sqrt(std/number_of_samples);
    return std;
}
/*
argv[]:
0: name of executable
1: input file name
2: output file order
*/
int main(int argc, char *argv[])
{
    #ifdef DEBUG
    printf("In main.\n");
    #endif
    if (argc == 3 || argc == 4)
        init(argv[1]);
    else if (argc == 17)
    {
        init_parameter(argv);
        data = loadcsv(argv[16]);
    }
    else
        fprintf(stderr, "wrong commond.\n");
    
    
    #ifdef DEBUG
        printf("k1 = %lf\nk2 = %lf\nk3 = %lf\nk4 = %lf\nk5 = %lf\nlength = %lf\n",
        k1, k2, k3, k4, k5, length);
        printf("nvar = %d\nx1 = %f\nx2 = %f\neps = %g\nh1 = %g\nhmin = %g\n",
    nvar, x1, x2, eps, h1, hmin);
        printf("\n");
    #endif
    if (argc == 4) data = loadcsv(argv[3]);
    clock_t start_odeint = clock();
    odeint(vstart, nvar, x1, x2, eps, h1, hmin, &nok, &nbad,
        derivs,
        rkqs);
    clock_t time_odeint= clock() - start_odeint;

    #ifdef DEBUG
        printf("odeint completed.\n");
    #endif
    if (!number_of_samples)
    {
        clock_t start_save = clock();
        save_results(argv[2]);
        clock_t time_save = clock() - start_save;
        fprintf(stdout, "%e\n", (double)(time_odeint + time_save)/CLOCKS_PER_SEC);
    }
    else
    {
        clock_t start_calculate_std = clock();
        double std = calculate_std();
        clock_t time_calculate_std = clock() - start_calculate_std;
        fprintf(stdout, "%e,%e\n", std, (double)(time_odeint + time_calculate_std)/CLOCKS_PER_SEC);
        free_matrix(data, 1, number_of_samples, 1, 6);
    }
    free_matrix(yp, 1, nvar, 1, kmax);
    free_vector(xp, 1, kmax);
    #ifdef DEBUG
        printf("Calculation complete!\n");
    #endif
    return 0;
}
