#include <cstring>
#include "../../includes/Utilities/Transpose.h"

void transpose(double *A, unsigned N)
{
    unsigned i,j;
    double B[N][N];
    #pragma omp parallel for private(j) shared(A,B)
    for(i=0;i<N;i++) 
        for(j=0;j<N;j++)
            B[j][i] =   A[i*N+j];

    memcpy(A,&B[0][0],N*N*sizeof(double));
    return ;
}

void transpose(double *A, double *B, unsigned N)
{
    unsigned i,j;
    #pragma omp parallel for private(j) shared(A,B)
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
            B[j*N+i] =   A[i*N+j];
            
    return ;
}
