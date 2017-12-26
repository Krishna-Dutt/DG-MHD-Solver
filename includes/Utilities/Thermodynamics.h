#ifndef THERMODYNAMICS_H
#define THERMODYNAMCIS_H

double Pressure(double Rho, double IE);
double SoundSpeed( double Rho, double Pr);
double InternalEnergy(double Rho, double Pr);

void KineticEnergy(double a, double* rho, double b, double* u, double c, double* v, unsigned index, unsigned size, double* ke);
void MomentumFluxPressure(double a, double* f, double b, double* u, double c, double* Pr, unsigned index, unsigned size, double* momflux_P);
void MomentumFlux(double a, double* f, double b, double* u, unsigned index, unsigned size, double* momflux);
void EnergyFlux(double a, double* E, double b, double* Pr, double c, double* u, unsigned index, unsigned size, double* eflux);
void IE(double a, double* rho, double b, double* Tp, double c, double* Pr, unsigned index, unsigned size, double* ie);
void Temperature(double a, double* ie, double b, double* rho, unsigned index, unsigned size, double* Tp);
void SoundSpeed(double a, double* rho, double b, double* Pr, unsigned index, unsigned size, double* c);
void Pressure(double a, double* rho, double b, double* ie, unsigned index, unsigned size, double* Pr);

#endif