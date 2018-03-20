#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/PolyEval.h"
#include "../../includes/Utilities/MaterialProperties.h"
 

double BoundaryDensity(double x, double y) {
    return 1.0;
} 

double BoundaryU(double x, double y){
    // if ( y > 1e-6) return 0.135;
  //return 0.0;
  if ( x < 1e-6 && y > 1e-6) return 0.135*1;//return (0.27/(0.025*0.025)) * y *( 0.05 - y);
  return 0.0;
}

double BoundaryV(double x, double y) {
    return 0.0;
}


double BoundaryPressure(double x, double y) {
     if( x < 1e-6) return 7.142857142857143e+01;
  return 70.78057142;
  //return 7.142857142857143e+01 + x * (69.70057142 - 7.142857142857143e+01 );
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
