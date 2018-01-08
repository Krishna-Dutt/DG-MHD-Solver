#ifndef MHDFUNCTIONS_H
#define MHDFUNCTIONS_H

void TotalPressureMHD(double a, double* p, double b, double* bx, double c, double* by, double d, double* bz, unsigned index, unsigned size, double* pt);
void MHDIE(double a, double *rho, double b, double *P, double c, double *bdotb, unsigned index, unsigned size, double *ie);
void Pressure_MHD(double a, double *rho, double b, double *ie, double c, double *bx, double d, double *by, double e, double *bz, unsigned index, unsigned size, double *P);
void MHDMomentumFluxPressure(double a, double *rhoU, double b, double *U, double c, double *Pt, double d, double *bx, unsigned index, unsigned size, double *Momx);
void MHDMomentumFlux(double a, double *rhoU, double b, double *V, double c, double *bx, double d, double *by, unsigned index, unsigned size, double *momflux);
void MHDEnergyFlux(double a, double *rhoE, double b, double *Pt, double c, double *U, double d, double *bx, double e, double *vdotb, unsigned index, unsigned size, double *magflux);
void MHDMagFieldFlux(double a, double *U, double b, double *by, double c, double *V, double d, double *bx, unsigned index, unsigned size, double *magflux);

void MHDMaxEigenValue(double a, double *rho, double b, double *P, double c, double *bx, double d, double *vdotb, unsigned index, unsigned size, double *C); 


#endif