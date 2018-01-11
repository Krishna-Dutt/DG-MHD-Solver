#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/PolyEval.h"
#include "../../includes/Utilities/VectorOperations.h"

void AdotB2d(double a, double *A1, double b, double *A2, double c, double *B1, double d, double *B2, unsigned index, unsigned size, double *C) {
    for(int i=0; i<size; i+= index) {
        C[i] = (a*A1[i] * c*B1[i]) + (b*A2[i] * d*B2[i]);
    }

    return ;
}
void AdotB3d(double a, double *A1, double b, double *A2, double c, double *A3, double d, double *B1, double e, double *B2, double f, double *B3, unsigned index, unsigned size, double *C) {
    for(int i=0; i<size; i+= index) {
        C[i] = (a*A1[i] * d*B1[i]) + (b*A2[i] * e*B2[i]) + (c*A3[i] * f*B3[i]);
    }

    return ;
}
