#include <stdio.h>
#include "nrutil.h"
#include "loadcsv.h"

const int MAX_SAMPLE_NUMBER = 10000;
const int NUMBER_OF_COLUMNS = 2;

extern int number_of_samples;
double **loadcsv(char* filename)
{
    //printf("csv file path is %s.\n", filename);
    FILE *file = fopen(filename, "r");
    //printf("file %s opened.\n", filename);
    char header_buffer[100];
    double **buffer = matrix(1, MAX_SAMPLE_NUMBER, 1, NUMBER_OF_COLUMNS);
    double **data;

    double time, absorbance;
    int counter = 0;
    fscanf(file, "%s\n", header_buffer);
    while (fscanf(file, "%lf,%lf\n", &time, &absorbance) != EOF)
    {
        counter++;
        buffer[counter][1] = time;
        buffer[counter][2] = absorbance;
        //printf("%lf,%lf\n", buffer[counter][0], buffer[counter][1]);
        //if(counter <10) printf("%lf,%lf\n", buffer[counter][0], buffer[counter][1]);
    }
    fclose(file);

    number_of_samples = counter;
    data = matrix(1, counter, 1, NUMBER_OF_COLUMNS);
    for (int i = 1; i <= counter; i++)
    {
        data[i][1] = buffer[i][1];
        data[i][2] = buffer[i][2];
        //if(i < 10) printf("%lf,%lf\n", data[i][0], data[i][1]);
    }
    free_matrix(buffer, 1, MAX_SAMPLE_NUMBER, 1, NUMBER_OF_COLUMNS);
    return data;
}