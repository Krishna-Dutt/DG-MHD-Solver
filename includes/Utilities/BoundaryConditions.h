#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

// Euler System 

double BoundaryDensity(double x, double y); 
double BoundaryU(double x, double y);
double BoundaryV(double x, double y);
double BoundaryPressure(double x, double y);
double BoundaryGamma(double x, double y);
double BoundaryEnergy(double x, double y);
double BoundarySpecificEnthalpy(double x, double y);
double BoundaryMachNo(double x, double y);
double BoundarySoundSpeed(double x, double y);

double BoundaryDensity(double x, double y, double t); 
double BoundaryU(double x, double y, double t);
double BoundaryV(double x, double y, double t);
double BoundaryPressure(double x, double y, double t);
double BoundaryGamma(double x, double y, double t);
double BoundaryEnergy(double x, double y, double t);
double BoundarySpecificEnthalpy(double x, double y, double t);
double BoundaryMachNo(double x, double y, double t);
double BoundarySoundSpeed(double x, double y, double t);

#endif
