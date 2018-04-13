#ifndef THERMODYNAMICS_H
#define THERMODYNAMICS_H

#define gamma (1.4)
#define R     287.15
#define Cp    (gamma*R/(gamma - 1.0))

// Navier-Stokes
#define meu 1e-3
#define PrandtlNo  0.72

// Sutherland's Coefficients
#define meu_ref 1.716e-5
#define T_ref   273.15
#define S_ref   110.4

// MHD
#define meu_perm 1.0

double Meu( double T);
double ThermalConductivity( double T);

#endif
