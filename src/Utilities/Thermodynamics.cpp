#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/PolyEval.h"
 
// Euler System 

double Pressure(double Rho, double IE) {
    return IE * (gamma - 1.0);
}

double SoundSpeed( double Rho, double P) {
    return sqrt(gamma * P/Rho);
}

double InternalEnergy(double Rho, double P) {
    return P/(gamma -1.0);
}