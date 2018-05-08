#ifndef TIMEINTEGRATORS_H
#define TIMEINTEGRATORS_H

void RK3step1(double t, double *dt, double a, double *A, double b, double *B, unsigned index, unsigned size, double *Z);
void RK3step2(double t, double *dt, double a, double *A, double b, double *B, double c, double *C, unsigned index, unsigned size, double *Z);
void RK3step3(double t, double *dt, double a, double *A, double b, double *B, double c, double *C, double d, double *D, unsigned index, unsigned size, double *Z);


#endif