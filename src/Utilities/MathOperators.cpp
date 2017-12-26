#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/MaterialProperties.h"
#include "../../includes/Utilities/MathOperators.h"

 
void Product(double a, double* x, double b, double* y, unsigned index, unsigned size, double * z) {
    for(int i=0; i< size; i+=index) {
      z[i] = a * x[i] * b * y[i];
    }

    return ;
}

void Addab(double a, double* x, double b, double* y, unsigned index, unsigned size, double * z) {
    for(int i=0; i< size; i+=index) {
      z[i] = (a * x[i]) + ( b * y[i]);
    }

    return ;
}

void Addabc(double a, double* x, double b, double* y, double c, double* z, unsigned index, unsigned size, double * Z) {
    for(int i=0; i< size; i+=index) {
      Z[i] = (a * x[i]) + ( b * y[i]) + (c * z[i]);
    }

    return ;
}

void Addabcd(double a, double* p, double b, double* q, double c, double* r, double d, double* s, unsigned index, unsigned size, double * t) {
    for(int i=0; i< size; i+=index) {
      t[i] = (a * p[i]) + ( b * q[i]) + (c * r[i]) + (d * s[i]);
    }

    return ;
}

void Subtract(double a, double* x, double b, double* y, unsigned index, unsigned size, double * z) {
    for(int i=0; i< size; i+=index) {
      z[i] = (a * x[i]) - (b * y[i]);
    }

    return ;
}

void Divide(double a, double* x, double b, double* y, unsigned index, unsigned size, double * z) {
    for(int i=0; i< size; i+=index) {
      if (b*y[i] != 0)
          z[i] = (a * x[i]) / (b * y[i]);
      else {
        std::cout << "Error :: Division by Zero !! \n";
      }
    }

    return ;
}
 
void ModulusAdd(double a, double* x, double b, double* y, unsigned index, unsigned size, double * z) {
    for(int i=0; i< size; i+=index) {
      z[i] = abs(a * x[i]) + abs( b * y[i]);
    }

    return ;
}


