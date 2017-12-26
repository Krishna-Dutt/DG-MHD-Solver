#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/PolyEval.h"
#include "../../includes/Utilities/ThermodynamicFunctions.h"
#include "../../includes/Utilities/MaterialProperties.h"

 
// Euler System 

double ReturnPressure(double Rho, double IE) {
    return IE * (gamma - 1.0);
}

double ReturnSoundSpeed( double Rho, double Pr) {
    return sqrt(gamma * Pr/Rho);
}

double ReturnInternalEnergy(double Rho, double Pr) {
    return Pr/(gamma -1.0);
}

void KineticEnergy(double a, double* rho, double b, double* u, double c, double* v, unsigned index, unsigned size, double* ke) {
  for(int i=0; i < size; i+=index) {
    ke[i] = 0.5 * (a*rho[i]) * (pow(b*u[i], 2) + pow(c*v[i],2));
  }
  return ;
}

void MomentumFluxPressure(double a, double* f, double b, double* u, double c, double* Pr, unsigned index, unsigned size, double* momflux_P) {
  for(int i=0; i < size; i+=index) {
    momflux_P[i] = (a*f[i]) * (b*u[i]) + (c*Pr[i]);
  }
  return ;
}
  

void MomentumFlux(double a, double* f, double b, double* u, unsigned index, unsigned size, double* momflux) {
  for(int i=0; i < size; i+=index) {
    momflux[i] = (a*f[i]) * (b*u[i]);
  }
  return ;
}


void EnergyFlux(double a, double* E, double b, double* Pr, double c, double* u, unsigned index, unsigned size, double* eflux) {
  for(int i=0; i < size; i+=index) {
    eflux[i] = ((a*E[i]) + (b*Pr[i])) * (c * u[i]) ;
  }
  return ;
}

void IE(double a, double* rho, double b, double* Tp, double c, double* Pr, unsigned index, unsigned size, double* ie) {
  for(int i=0; i < size; i+=index) {
    ie[i] = Pr[i]/(gamma - 1.0) ;
  }
  return ;
}
  
void Temperature(double a, double* ie, double b, double* rho, unsigned index, unsigned size, double* Tp) {
  for(int i=0; i < size; i+=index) {
    Tp[i] = ie[i]*( gamma - 1.0 )/(rho[i] * R) ;
  }
  return ;
}

void SoundSpeed(double a, double* rho, double b, double* Pr, unsigned index, unsigned size, double* c) {
  for(int i=0; i < size; i+=index) {
    c[i] = sqrt(gamma * Pr[i]/ rho[i]) ;
  }
  return ;
}

void Pressure(double a, double* rho, double b, double* ie, unsigned index, unsigned size, double* Pr) {
    for(int i=0; i < size; i+=index) {
    Pr[i] = ie[i] * (gamma - 1.0) ;
  }
  return ;
}

