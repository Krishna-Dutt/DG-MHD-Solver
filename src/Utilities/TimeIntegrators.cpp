#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/MaterialProperties.h"
#include "../../includes/Utilities/TimeIntegrators.h"

void RK3step1(double t, double *dt, double a, double *A, double b, double *B, unsigned index, unsigned size, double *Z){
    
    #pragma omp parallel for
    for(int i=0; i< size; i+=index) {
      Z[i] = t * dt[i] * ( a * A[i] ) + b * B[i]  ;
    }

    return ;
}

void RK3step2(double t, double *dt, double a, double *A, double b, double *B, double c, double *C, unsigned index, unsigned size, double *Z){
    
    #pragma omp parallel for
    for(int i=0; i< size; i+=index) {
      Z[i] = t * dt[i] * ( a * A[i] + b * B[i] ) + c * C[i] ;
    }

    return ;
}
void RK3step3(double t, double *dt, double a, double *A, double b, double *B, double c, double *C, double d, double *D, unsigned index, unsigned size, double *Z) {

    #pragma omp parallel for
    for(int i=0; i< size; i+=index) {
      Z[i] = t * dt[i] * ( a * A[i] + b * B[i] + c * C[i] ) + d * D[i] ;
    }

    return ;
}
