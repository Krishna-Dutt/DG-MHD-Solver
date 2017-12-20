#ifndef LOBATTOINTEGRATION_H
#define LOBATTOINTEGRATION_H

#include <functional>

using namespace std;

double lobattoIntegration(double start, double end, unsigned N, function<double(double)> f);
double lobattoIntegration(double start, double end, int N, double **f_value);
double lobattoIntegration(double start, double end, int N, unsigned Index, double *f_value);


#endif
