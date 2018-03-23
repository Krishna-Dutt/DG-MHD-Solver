#ifndef THERMODYNAMICFUNCTIONS_H
#define THERMODYNAMICFUNCTIONS_H

double ReturnPressure(double Rho, double IE);
double ReturnSoundSpeed( double Rho, double Pr);
double ReturnInternalEnergy(double Rho, double Pr);

void KineticEnergy(double a, double* rho, double b, double* u, double c, double* v, unsigned index, unsigned size, double* ke);
void KineticEnergy3d(double a, double* rho, double b, double* u, double c, double* v, double d, double* w, unsigned index, unsigned size, double* ke);
void MomentumFluxPressure(double a, double* f, double b, double* u, double c, double* Pr, unsigned index, unsigned size, double* momflux_P);
void MomentumFlux(double a, double* f, double b, double* u, unsigned index, unsigned size, double* momflux);
void EnergyFlux(double a, double* E, double b, double* Pr, double c, double* u, unsigned index, unsigned size, double* eflux);
void IE(double a, double* rho, double b, double* Tp, double c, double* Pr, unsigned index, unsigned size, double* ie);
void Temperature(double a, double* ie, double b, double* rho, unsigned index, unsigned size, double* Tp);
void SoundSpeed(double a, double* rho, double b, double* Pr, unsigned index, unsigned size, double* c);
void Pressure(double a, double* rho, double b, double* ie, unsigned index, unsigned size, double* Pr);
void ThermoEntropy(double a, double* Pr, double b, double* rho, unsigned index, unsigned size, double* en);
void EntropyVar(double a, double* rho, double b, double* Pr, double c, double* BdotB, unsigned index, unsigned size, double* en);

void ViscousStress(double a,  double* x, double b, double* y, unsigned index, unsigned size, double* z);
void EnergyViscous(double a, double* u, double b, double* uTau, double c, double* v, double d, double* vTau, unsigned index, unsigned size, double* Evisc);
void HeatFlux(double a, double *T, double b, double *gradT, unsigned index, unsigned size, double *z);
#endif