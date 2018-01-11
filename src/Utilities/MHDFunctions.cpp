#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/PolyEval.h"
#include "../../includes/Utilities/MHDFunctions.h"
#include "../../includes/Utilities/MaterialProperties.h"


void TotalPressureMHD(double a, double* p, double b, double* bx, double c, double* by, double d, double* bz, unsigned index, unsigned size, double* pt){
    for(int i=0; i< size; i+=index) {
        pt[i] = a*p[i] + 0.5 * (pow(b*bx[i], 2) + pow(c*by[i], 2) + pow(d*bz[i], 2));
    }

    return ;
}

void MHDIE(double a, double *rho, double b, double *P, double c, double *bdotb, unsigned index, unsigned size, double *ie){
    for(int i=0; i<size; i+= index) {
        ie[i] = a*P[i]/(gamma -1.0)  + 0.5 * bdotb[i];
    }

    return ;
}

void Pressure_MHD(double a, double *rho, double b, double *ie, double c, double *bx, double d, double *by, double e, double *bz, unsigned index, unsigned size, double *P){
    for(int i=0; i<size; i+=index) {
        P[i] = (gamma -1.0) * ( b*ie[i] - 0.5 * (pow(c*bx[i],2) + pow(d*by[i],2) + pow(e*bz[i],2) )) ;
    }

    return ;
}

void MHDMomentumFluxPressure(double a, double *rhoU, double b, double *U, double c, double *Pt, double d, double *bx, unsigned index, unsigned size, double *Momx){
    for(int i=0; i<size; i+= index) {
        Momx[i] = a*rhoU[i] * b*U[i] + c*Pt[i] - pow(d*bx[i],2);
    }

    return ;
}

void MHDMomentumFlux(double a, double *rhoU, double b, double *V, double c, double *bx, double d, double *by, unsigned index, unsigned size, double *momflux){
    for(int i=0; i<size; i+=index) {
        momflux[i] = a*rhoU[i] * b*V[i]  - c*bx[i] * d*by[i];
    }

    return ;
}

void MHDEnergyFlux(double a, double *rhoE, double b, double *Pt, double c, double *U, double d, double *bx, double e, double *vdotb, unsigned index, unsigned size, double *magflux){
    for(int i=0; i<size; i+= index) {
        magflux[i] = ( a*rhoE[i] + b*Pt[i]) * c*U[i] - d*bx[i] * e*vdotb[i];
    }

    return ;
}

void MHDMagFieldFlux(double a, double *U, double b, double *by, double c, double *V, double d, double *bx, unsigned index, unsigned size, double *magflux){
    for(int i=0; i<size; i+=index) {
        magflux[i] = a*U[i] * b*by[i] - c*V[i] * d*bx[i];
    }
    
    return ;
}

void MHDMaxEigenValue(double a, double *rho, double b, double *P, double c, double *bx, double d, double *bdotb, unsigned index, unsigned size, double *C){
    double temp_a = 0 , temp_b = 0, temp_bx = 0;
    for(int i=0; i<size; i+=index){
        temp_a = gamma * P[i]/rho[i];
        temp_bx = bx[i]/sqrt(rho[i]);
        temp_b = bdotb[i]/rho[i];
        //std::cout << temp_a << " : " << temp_b << " : " << temp_bx << endl;
        C[i] = sqrt(0.5 * (temp_a + temp_b + sqrt(pow(temp_a + temp_b,2) -4.0*temp_a*temp_bx*temp_bx)));
    }

    return ;
}
