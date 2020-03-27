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
    clock_t start_odeint, end_odeint, start_save, end_save;
    double time_odeint, time_save;
    #ifdef DEBUG
    printf("In main.\n");
    #endif
    printf("argc == %d\n", argc);
    if (argc == 3) // no experiment data to be compared
    {
        init(argv[1]);
        #ifdef DEBUG
            printf("k1 = %lf\nk2 = %lf\nk3 = %lf\nk4 = %lf\nk5 = %lf\nlength = %lf\n",
            k1, k2, k3, k4, k5, length);
            printf("nvar = %d\nx1 = %f\nx2 = %f\neps = %g\nh1 = %g\nhmin = %g\n",
        nvar, x1, x2, eps, h1, hmin);
            printf("\n");
        #endif
        start_odeint = clock();
        odeint(vstart, nvar, x1, x2, eps, h1, hmin, &nok, &nbad,
            derivs,
            rkqs);
        end_odeint = clock();
        time_odeint = (double)(end_odeint - start_odeint)/CLOCKS_PER_SEC;

        #ifdef DEBUG
            printf("odeint completed.\n");
        #endif
        start_save = clock();
        save_results(argv[2]);
        end_save = clock();
        time_save = (double)(end_save-start_save)/CLOCKS_PER_SEC;
        free_space();
        char path[]="nok_and_nbad.txt\0";
        FILE *file = fopen(path, "w");
        fprintf(file, "%d\n%d\n%lf\n%lf", nok, nbad, time_odeint, time_save);
        fclose(file);
        #ifdef DEBUG
            printf("Calculation complete!\n");
        #endif
    }

    if (argc == 4) // need to compare with experiment data and consider timestep fitting during simulation
    {
        int number_of_samples = 0;
        data = loadcsv(argv[3], &number_of_samples);
        printf("number of samples is %d.\n", number_of_samples);
        printf("address of data in main is %d.\n", data);
        for (int i = 1; i <= number_of_samples; i++)
            printf("%lf,%lf\n", data[i][1], data[i][2]);
        
        free_matrix(data, 1, number_of_samples, 1, 2);
        
    }
    return 0;
}
