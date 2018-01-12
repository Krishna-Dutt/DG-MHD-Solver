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

// MHD
 double BoundaryMHDEnergy(double x, double y);
 double BoundaryBX(double x, double y);
 double BoundaryBY(double x, double y);
 double BoundaryBZ(double x, double y);
 double BoundaryW(double x, double y);

 double BoundaryMHDEnergy(double x, double y, double t);
 double BoundaryBX(double x, double y, double t);
 double BoundaryBY(double x, double y, double t);
 double BoundaryBZ(double x, double y, double t);
 double BoundaryW(double x, double y, double t);
  

#endif
