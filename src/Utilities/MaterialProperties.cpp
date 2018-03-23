#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/MaterialProperties.h"

double Meu( double T) {

    return meu_ref * pow(T/T_ref, 1.5) * (T_ref + S_ref)/(T + S_ref);
}

double ThermalConductivity(double T) {
    
    return Cp * Meu(T)/PrandtlNo;
}