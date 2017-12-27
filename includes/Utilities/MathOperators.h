#ifndef MATHOPERATORS_H
#define MATHOPERATORS_H

void Product(double a, double* x, double b, double* y, unsigned index, unsigned size, double * z);
void Addab(double a, double* x, double b, double* y, unsigned index, unsigned size, double * z);
void Addabc(double a, double* x, double b, double* y, double c, double* z, unsigned index, unsigned size, double * Z);
void Addabcd(double a, double* p, double b, double* q, double c, double* r, double d, double* s, unsigned index, unsigned size, double * t);
void Subtract(double a, double* x, double b, double* y, unsigned index, unsigned size, double * z);
void Divide(double a, double* x, double b, double* y, unsigned index, unsigned size, double * z);
void ModulusAdd(double a, double* x, double b, double* y, unsigned index, unsigned size, double * z);

void Copy(double a, double* x, unsigned index, unsigned size, double* z);
#endif