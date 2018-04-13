#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/PolyEval.h"
#include "../../includes/Utilities/MaterialProperties.h"
 

double BoundaryDensity(double x, double y) {
    if ( y <= -tan(30.8*M_PI/180.0)*(x -1.0)) return 2.445e-3;
  return 2.862e-3;
} 

double BoundaryU(double x, double y){
  if ( y <= -tan(30.8*M_PI/180.0)*(x -1.0) ) return 731.756;
  return 705.355*cos(3.813*M_PI/180.0);
}

double BoundaryV(double x, double y) {
    if ( y <= -tan(30.8*M_PI/180.0)*(x -1.0) ) return 0.0;
  return -705.355*sin(3.813*M_PI/180.0);
}


double BoundaryPressure(double x, double y) {
     if ( y <= -tan(30.8*M_PI/180.0)*(x -1.0)) return 2.445e-3*R*288.15;
  return 2.862e-3*R*307.025;
}

double BoundaryTemperature(double x, double y) {
    if ( y <= -tan(30.8*M_PI/180.0)*(x -1.0) ) return 288.15;
  return 307.025;
}

double BoundaryGamma(double x, double y) {
    return gamma;
}
    
double BoundaryEnergy(double x, double y) {
    return BoundaryPressure(x,y)/(BoundaryGamma(x,y) -1.0) + 0.5 * BoundaryDensity(x,y) *(pow(BoundaryU(x,y),2) +pow(BoundaryV(x,y),2));
}

double BoundarySpecificEnthalpy(double x, double y) {
    return  BoundaryGamma(x,y)*BoundaryPressure(x,y)/((BoundaryGamma(x,y) -1.0)*BoundaryDensity(x,y)) + 0.5 *(pow(BoundaryU(x,y),2) +pow(BoundaryV(x,y),2));
}

double BoundarySoundSpeed(double x, double y) {
    return sqrt( BoundaryGamma(x,y) * BoundaryPressure(x,y) / BoundaryDensity(x,y) );
}

double BoundaryMachNo(double x, double y) {
    return sqrt(pow(BoundaryU(x,y),2) +pow(BoundaryV(x,y),2))/BoundarySoundSpeed(x,y);
}

double BoundaryDensity(double x, double y, double t) {
     if( x < 0.5) return 1.0;
     return 0.125;
} 

double BoundaryU(double x, double y, double t){
    return 0.0;
}

double BoundaryV(double x, double y, double t) {
    return 0.0;
}

double BoundaryPressure(double x, double y, double t) {
     return 1.0;
}

double BoundaryTemperature(double x, double y, double t) {
     return 1.0;
}

double BoundaryGamma(double x, double y, double t) {
    return gamma;
}
    
double BoundaryEnergy(double x, double y, double t) {
    return BoundaryPressure(x,y,t)/(BoundaryGamma(x,y,t) -1.0) + 0.5 * BoundaryDensity(x,y,t) *(pow(BoundaryU(x,y,t),2) +pow(BoundaryV(x,y,t),2));
}

double BoundarySpecificEnthalpy(double x, double y, double t) {
    return  BoundaryGamma(x,y,t)*BoundaryPressure(x,y,t)/((BoundaryGamma(x,y,t) -1.0)*BoundaryDensity(x,y,t)) + 0.5 *(pow(BoundaryU(x,y,t),2) +pow(BoundaryV(x,y,t),2));
}

double BoundarySoundSpeed(double x, double y, double t) {
    return sqrt( BoundaryGamma(x,y,t) * BoundaryPressure(x,y,t) / BoundaryDensity(x,y,t) );
}

double BoundaryMachNo(double x, double y, double t) {
    return sqrt(pow(BoundaryU(x,y,t),2) +pow(BoundaryV(x,y,t),2))/BoundarySoundSpeed(x,y,t);
}


// MHD 

double BoundaryBX(double x, double y){
    return 1.0;
}

double BoundaryBY(double x, double y) {
    return 1.0;
}

double BoundaryBZ(double x, double y) {
    return 1.0;
}

double BoundaryW(double x, double y) {
    return 1.0;
}


double BoundaryW(double x, double y, double t) {
    return 1.0;
}

double BoundaryBX(double x, double y, double t) {
    return 1.0;
}

double BoundaryBY(double x, double y, double t) {
    return 1.0;
}

double BoundaryBZ(double x, double y, double t) {
    return 1.0;
}

double BoundaryMHDEnergy(double x, double y) {
    return BoundaryPressure(x,y)/(BoundaryGamma(x,y) -1.0) + 0.5 * BoundaryDensity(x,y) *(pow(BoundaryU(x,y),2) +pow(BoundaryV(x,y),2)) + 0.5*(pow(BoundaryBX(x,y),2) +pow(BoundaryBY(x,y),2) +pow(BoundaryBZ(x,y),2));
}

double BoundaryMHDEnergy(double x, double y, double t) {
    return BoundaryPressure(x,y,t)/(BoundaryGamma(x,y,t) -1.0) + 0.5 * BoundaryDensity(x,y,t) *(pow(BoundaryU(x,y,t),2) +pow(BoundaryV(x,y,t),2)) + 0.5*(pow(BoundaryBX(x,y,t),2) +pow(BoundaryBY(x,y,t),2) +pow(BoundaryBZ(x,y,t),2));
}
