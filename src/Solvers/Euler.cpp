#include "../../includes/Solvers/EulerSolver.h"
#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/MaterialProperties.h"
#include "../../includes/Utilities/MathOperators.h"
#include "../../includes/Utilities/ThermodynamicFunctions.h"


EulerSolver::EulerSolver(int _ne_x, int _ne_y, int _N) {
    ne_x = _ne_x;
    ne_y = _ne_y;
    N = _N;
    time = 0.0;
}

EulerSolver::~EulerSolver() {
  delete field;
}

void EulerSolver::setDomain(double _x1, double _y1, double _x2, double _y2) {
    x1 = _x1;
    y1 = _y1;
    x2 = _x2;
    y2 = _y2;
    field = new DG_Field_2d(ne_x, ne_y, N, x1, y1, x2, y2);
    field->setSystem("EULER");

   return ;
}

void EulerSolver::setPrimitiveVariables(){
  D  = field->addVariable_withBounary("q");
  Vx = field->addVariable_withBounary("u");
  Vy = field->addVariable_withBounary("v");
  P  = field->addVariable_withBounary("P");
  T  = field->addVariable_withBounary("T");
  return ;
}

/*void EulerSolver::setGradientPrimitiveVariables(){
  field->addVariable_withoutBounary("dqdx");
  field->addVariable_withoutBounary("dudx");
  field->addVariable_withoutBounary("dvdx");
  field->addVariable_withoutBounary("dPdx");
  field->addVariable_withoutBounary("dTdx");

  field->addVariable_withoutBounary("dqdy");
  field->addVariable_withoutBounary("dudy");
  field->addVariable_withoutBounary("dvdy");
  field->addVariable_withoutBounary("dPdy");
  field->addVariable_withoutBounary("dTdy");

  return ;
}*/


void EulerSolver::setConservativeVariables(){
  DVx = field->addVariable_withBounary("qu");
  DVy = field->addVariable_withBounary("qv");
  DE  = field->addVariable_withBounary("qE");
  De  = field->addVariable_withBounary("qe"); // Internal Energy, added just for ease of manipulating Energy, Remove later if not required !!
  KE  = field->addVariable_withBounary("KE"); // Kinetic Energy , " "

  field->addConservativeVariables(D);
  field->addConservativeVariables(DVx);
  field->addConservativeVariables(DVy);
  field->addConservativeVariables(DE);
  return ;
}

/*void EulerSolver::setMaterialPropertyVariables(){
  field->addVariable_withoutBounary("meu");
  return ;
}*/

void EulerSolver::setInitialDensity(function<double(double, double)> Rho) {
    field->initializeVariable(D, Rho);
    return ;
}

void EulerSolver::setInitialPressure(function<double(double, double)> Pr) {
    field->initializeVariable(P, Pr);
    return ;
}

void EulerSolver::setInitialTemperature(function<double(double, double)> Tp) {
    field->initializeVariable(T, Tp);
    return ;
}

void EulerSolver::setInitialVelocity(function<double(double, double)> U, function<double(double, double)> V) {
    field->initializeVariable(Vx, U);
    field->initializeVariable(Vy, V);
    return ;
}

void EulerSolver::setBoundaryCondtions(string type1, string type2, string type3, string type4) {
      // Setting BC as Outflow type to test Methods
    field->setBottomBoundary(type1);
    field->setRightBoundary(type2);
    field->setTopBoundary(type3);
    field->setLeftBoundary(type4);

    return ;
}

void EulerSolver::setSolver(double _CFL, double _time, int _no_of_time_steps) {
   CFL = _CFL;
   time = _time;
   no_of_time_steps = _no_of_time_steps;
   Dx   = field->addVariable_CellCentered();
   field->FindMindx(Dx);
   Dt   = field->addVariable_CellCentered();
   UMax = field->addVariable_CellCentered();
   field->ResetVariables_CellCentered(UMax, 0.0);

    return ;
}

/*double Subtract2(double x, double y) {
  return (x ) ;
}*/


/*double PressureE(double D, double e) {
  return D*e*(1.4 -1.0) ;
}*/

/*double NormalViscousStress(double meu, double x, double y) {
  return meu * ((4.0/3.0)*x - (2.0/3.0)*y) ;
}*/

/*double TangentialViscousStress(double meu, double x, double y) {
  return meu * ( x + y) ;
}*/

/*double EnergyViscous(double u, double uTau, double v, double vTau) {
  return (u * uTau + v * vTau) ;
}*/

/*double ArtificialViscosity( double x , double y) {
  double Beta = 1.5;
  double B1 = 6.0;
  double B2 = 20.5;

 // return (Beta - 1.0)*(B1/B2)* x * (1.0/500) ;
  return 1e-4;//2e-5;
}*/

// Functions for cell centered variables
void Maximum( double a, double* x, double b, double* y, unsigned size_DV, unsigned size_CV, double* z, unsigned node) {
  int j =-1 , count = 0; 
  for(int i=0; i < size_DV; ++i) {
    if(count % node == 0) {
      j++;
      count = 0;
      z[j] = a*x[i];
    }
    z[j] = max( max(a*x[i], b*y[i]), z[j]);
    count++;
  }

  return ;
}

void EulerSolver::setXMomentum() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vx, Product, DVx);
  return ;
}

void EulerSolver::setYMomentum() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vy, Product, DVy);
  return ;
}

void EulerSolver::setEnergy() {
  field->setFunctionsForVariables(1.0, KE, 1.0, De, Addab, DE);
  return ;
}

void EulerSolver::setInternalEnergyfromPrimitive() {
  field->setFunctionsForVariables(1.0, D, 1.0, T, 1.0, P, IE, De);
  return ;
}

void EulerSolver::setInternalEnergy() {
  field->setFunctionsForVariables(1.0, DE, 1.0, KE, Subtract, De);
  return ;
}


void EulerSolver::setKineticEnergy() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vx, 1.0, Vy, KineticEnergy, KE);
  return ;
}

void EulerSolver::updateVelocity() {
  field->setFunctionsForVariables(1.0, DVx, 1.0, D, Divide, Vx);
  field->setFunctionsForVariables(1.0 ,DVy, 1.0, D, Divide, Vy);
  return ;
}

void EulerSolver::updateTemperature() {
  field->setFunctionsForVariables(1.0, De, 1.0, D, Temperature, T);
  return ;
}

void EulerSolver::updatePressure() {
  field->setFunctionsForVariables(1.0, D, 1.0, De, Pressure, P);
  return ;
}

void EulerSolver::updatePrimitiveVariables() {
  updateVelocity();
  setKineticEnergy();
  setInternalEnergy();
  updateTemperature();
  updatePressure();
  return ;
}

void EulerSolver::updateConservativeVariables() {
  setXMomentum();
  setYMomentum();
  setKineticEnergy();
  setInternalEnergyfromPrimitive();
  setEnergy();
  return ;
}

void EulerSolver::setInviscidFlux() {
  //cout << " Calling setInviscidFlux " << endl;
  DVxVx_plus_P = field->addVariable_withBounary("quu_plus_P");
  DVxVy        = field->addVariable_withBounary("quv");
  DVyVy_plus_P = field->addVariable_withBounary("qvv_plus_P");
  DE_plus_P_Vx = field->addVariable_withBounary("qE_plus_P_u");
  DE_plus_P_Vy = field->addVariable_withBounary("qE_plus_P_v");

  updateInviscidFlux();
  return ;
}

void EulerSolver::updateInviscidFlux() {
  field->setFunctionsForVariables(1.0, DVx, 1.0, Vx, 1.0, P, MomentumFluxPressure, DVxVx_plus_P);
  field->setFunctionsForVariables(1.0, DVy, 1.0, Vy, 1.0, P, MomentumFluxPressure, DVyVy_plus_P);
  field->setFunctionsForVariables(1.0, DVx, 1.0, Vy, MomentumFlux, DVxVy);
  field->setFunctionsForVariables(1.0, DE, 1.0, P, 1.0, Vx, EnergyFlux, DE_plus_P_Vx);
  field->setFunctionsForVariables(1.0, DE, 1.0, P, 1.0, Vy, EnergyFlux, DE_plus_P_Vy);
  
  return ;
}

/*void EulerSolver::setViscousFlux() {
  field->addVariable_withBounary("Tauxx");
  field->addVariable_withBounary("Tauxy");
  field->addVariable_withBounary("Tauyy");
  field->addVariable_withBounary("Eviscousx");
  field->addVariable_withBounary("Eviscousy");

  updateInviscidFlux();
  return ;
}*/

/*void EulerSolver::updateViscousFlux() {
  field->setFunctionsForVariables("meu", "dudx", "dvdy", NormalViscousStress, "Tauxx");
  field->setFunctionsForVariables("meu", "dvdy", "dudx", NormalViscousStress, "Tauyy");
  field->setFunctionsForVariables("meu", "dudy", "dvdx", TangentialViscousStress, "Tauxy");
  field->setFunctionsForVariables(Vx, "Tauxx", Vy, "Tauxy", EnergyViscous, "Eviscousx");
  field->setFunctionsForVariables(Vx, "Tauxy", Vy, "Tauyy", EnergyViscous, "Eviscousy");
  return ;
}*/

/*void EulerSolver::updatePrimitiveGradient() {
  field->delByDelX(D, "dqdx", "central");
  field->delByDelX(Vx, "dudx", "central");
  field->delByDelX(Vy, "dvdx", "central");
  field->delByDelX(P, "dPdx", "central");
  field->delByDelX(T, "dTdx", "central");

  field->delByDelY(D, "dqdy", "central");
  field->delByDelY(Vx, "dudy", "central");
  field->delByDelY(Vy, "dvdy", "central");
  field->delByDelY(P, "dPdy", "central");
  field->delByDelY(T, "dTdy", "central");
}*/

void EulerSolver::setAuxillaryVariables() {
 // cout << " Calling setAuxillaryVariables " << endl;
  K1D   = field->addVariable_withoutBounary();
  K1DVx = field->addVariable_withoutBounary();
  K1DVy = field->addVariable_withoutBounary();
  K1DE  = field->addVariable_withoutBounary();
  K2D   = field->addVariable_withoutBounary();
  K2DVx = field->addVariable_withoutBounary();
  K2DVy = field->addVariable_withoutBounary();
  K2DE  = field->addVariable_withoutBounary();
  K3D   = field->addVariable_withoutBounary();
  K3DVx = field->addVariable_withoutBounary();
  K3DVy = field->addVariable_withoutBounary();
  K3DE  = field->addVariable_withoutBounary();
  dbydx = field->addVariable_withoutBounary();
  dbydy = field->addVariable_withoutBounary();

  DAnalytical  = field->addVariable_withBounary("qAnalytic");
  VxAnalytical = field->addVariable_withBounary("uAnalytic");
  
  ZERO = field->addVariable_withoutBounary();
  
  return ;
}

void EulerSolver::setEigenValues() {
 // cout << "Calling setEigenValues " << endl;
  C         = field->addVariable_withBounary("c");
  Vx_plus_C = field->addVariable_withBounary("u_plus_c");
  Vy_plus_C = field->addVariable_withBounary("v_plus_c");// Recheck formulation of eigen value !!

  updateEigenValues();
  return ;
}

void EulerSolver::updateEigenValues() {
 // cout << " Calling updateEigenValues " << endl;
  field->setFunctionsForVariables(1.0, D, 1.0, P, SoundSpeed, C);
  field->setFunctionsForVariables(1.0, Vx, 1.0, C, ModulusAdd, Vx_plus_C);
  field->setFunctionsForVariables(1.0, Vy, 1.0, C, ModulusAdd, Vy_plus_C);
  return ;
}

void EulerSolver::RK_Step1(int Var, int FluxX, int FluxY, int K) {
  field->delByDelX(FluxX, dbydx, Var, "rusanov", Vx_plus_C);
  field->delByDelY(FluxY, dbydy, Var, "rusanov", Vy_plus_C);

  field->setFunctionsForVariables(-1.0, dbydx, -1.0, dbydy, Addab, K);
  
  field->setFunctionsForVariables(0.5*dt, K, 1.0, Var, Addab, Var);

  return;
}

void EulerSolver::RK_Step2(int Var, int FluxX, int FluxY, int K1, int K2) {
  field->delByDelX(FluxX, dbydx, Var, "rusanov", Vx_plus_C);
  field->delByDelY(FluxY, dbydy, Var, "rusanov", Vy_plus_C);

  field->setFunctionsForVariables(-1.0, dbydx, -1.0, dbydy, Addab, K2);
  
  field->setFunctionsForVariables(-1.5*dt, K1, 2.0*dt, K2, 1.0, Var, Addabc, Var);

  return;
}

void EulerSolver::RK_Step3(int Var, int FluxX, int FluxY, int K1, int K2, int K3) {
  field->delByDelX(FluxX, dbydx, Var, "rusanov", Vx_plus_C);
  field->delByDelY(FluxY, dbydy, Var, "rusanov", Vy_plus_C);

  field->setFunctionsForVariables(-1.0, dbydx, -1.0, dbydy, Addab, K3);
  
  field->setFunctionsForVariables((7.0/6.0)*dt, K1, -(4.0/3.0)*dt, K2, (1.0/6.0)*dt, K3, 1.0, Var, Addabcd, Var);
  return;
}

void EulerSolver::setTimeStep() {
  field->setFunctionsForCellCenterVariablesfromDomainVariables(1.0, Vx_plus_C, 1.0, Vy_plus_C, Maximum, UMax);
  field->FindTimestep(Dt, Dx, UMax, CFL);
  dt = field->FindMindt(Dt);

  return ;
}


void EulerSolver::solve() {
  // Requires all Primitive and Conservative Variables to be setup and initialised.
  double t = 0.0;
  int count = 0;
  
 
  // Till Now all variables have to be initialised !!
  // For loop to march in time !!
  while(t <= time) {
    // First Step of RK3
    
    updateInviscidFlux();
    updateEigenValues();

    if ( count%no_of_time_steps == 0) {
      setTimeStep();
      cout << "Time Step : " << dt << " , Time : " << t << "\n"; 
      count = 0;
    }
    // Mass
    RK_Step1(D, DVx, DVy, K1D);
    // X Momentum
    RK_Step1(DVx,DVxVx_plus_P,DVxVy,K1DVx);
    // Y Momentum
    RK_Step1(DVy,DVxVy,DVyVy_plus_P,K1DVy);
    // Energy
    RK_Step1(DE,DE_plus_P_Vx,DE_plus_P_Vy,K1DE);
    
    RunShockDetector();
    RunLimiter();
    RunPositivityLimiter();

    field->updateBoundary(t);
    updatePrimitiveVariables();
   
    
    // Second Step of RK3
    updateInviscidFlux();
    updateEigenValues();

    // Mass
    RK_Step2(D, DVx, DVy, K1D, K2D);
    // X Momentum
    RK_Step2(DVx,DVxVx_plus_P,DVxVy,K1DVx, K2DVx);
    // Y Momentum
    RK_Step2(DVy,DVxVy,DVyVy_plus_P,K1DVy, K2DVy);
    // Energy
    RK_Step2(DE,DE_plus_P_Vx,DE_plus_P_Vy,K1DE, K2DE);
    
    RunShockDetector();
    RunLimiter();
    RunPositivityLimiter(); 
    
    field->updateBoundary(t);
    updatePrimitiveVariables();

   // Third (Final) Step of RK3
    updateInviscidFlux();
    updateEigenValues();

    // Mass
    RK_Step3(D, DVx, DVy, K1D, K2D, K3D);
    // X Momentum
    RK_Step3(DVx,DVxVx_plus_P,DVxVy,K1DVx, K2DVx, K3DVx);
    // Y Momentum
    RK_Step3(DVy,DVxVy,DVyVy_plus_P,K1DVy, K2DVy, K3DVy);
    // Energy
    RK_Step3(DE,DE_plus_P_Vx,DE_plus_P_Vy,K1DE, K2DE, K3DE);
 
    
   RunShockDetector();
   RunLimiter();
   RunPositivityLimiter(); 
   
    field->updateBoundary(t);
    updatePrimitiveVariables();

   t += dt; 
   count += 1;       
    }

  return ;
}

void EulerSolver::plot(string filename) {
    field->writeVTK(filename);
    return ;
}

void EulerSolver::SetShockDetector(string _ShockDetector) {
  ShockDetector = _ShockDetector;
  SetShockDetectorVariables();
  return ;
}

void EulerSolver::SetShockDetectorVariables() {
  
  if (ShockDetector == "KXRCF") {
    CellMarker  = field->addVariable_CellCentered();
    CellMarkerG = field->addVariable_withBounary("cellmarker");

    field->scal(0.0, CellMarkerG);
    
  }

  return ;
}

void EulerSolver::RunShockDetector() {
  if (ShockDetector == "KXRCF") {
    Run_KXRCF();
 }

  return ;
}

void EulerSolver::Run_KXRCF() {
 
  field->ResetVariables_CellCentered(CellMarker, 1.5);
  field->ResetMap_OutFlow();

  //field->updateOutFlowBoundary(Vx, Vy);
  //field->updateCellMarker(D, CellMarker);

  return ;
}



void EulerSolver::SetLimiter(string _Limiter) {
  Limiter = _Limiter;
  SetLimiterVariables();
  return ;
}

void EulerSolver::SetLimiterVariables() {
  if (Limiter == "LiliaMoment") {
    /*CellMarker  = field->addVariable_CellCentered();
    field->ResetVariables_CellCentered(CellMarker, 1.5);
    CellMarkerG = field->addVariable_withBounary("cellmarker");
    field->scal(0.0, CellMarkerG);*/

    Moment    = field->addVariable_withBounary("moment");
    field->scal(0.0, Moment);
    ModMoment = field->addVariable_withBounary("modmoment");
    field->scal(0.0, ModMoment);
    field->setVanderMandMatrix();
  }
  else if (Limiter == "CharacteristicLimiter") {
    field->setEigenMatrices(4);
    field->setVanderMandMatrix();
    uMoment = field->addVariable_withoutBounary();
    vMoment = field->addVariable_withoutBounary();
    qMoment = field->addVariable_withoutBounary();
    PMoment = field->addVariable_withoutBounary();
    HMoment = field->addVariable_withoutBounary();
    field->scal(0.0, uMoment);
    field->scal(0.0, vMoment);
    field->scal(0.0, qMoment);
    field->scal(0.0, PMoment);
    field->scal(0.0, HMoment);

    Char1 = field->addVariable_withoutBounary();
    Char2 = field->addVariable_withoutBounary();
    Char3 = field->addVariable_withoutBounary();
    Char4 = field->addVariable_withoutBounary();
    field->scal(0.0, Char1);
    field->scal(0.0, Char2);
    field->scal(0.0, Char3);
    field->scal(0.0, Char4);

    dPdx = field->addVariable_withBounary("dPdx");
    dPdy = field->addVariable_withBounary("dPdy");
    dPdxMoment = field->addVariable_withoutBounary();
    dPdyMoment = field->addVariable_withoutBounary();
    field->scal(0.0, dPdxMoment);
    field->scal(0.0, dPdyMoment);


    AuxillaryVariables["V"].resize(0);
    AuxillaryVariables["V"].push_back("qMoment");
    AuxillaryVariables["V"].push_back("uMoment");
    AuxillaryVariables["V"].push_back("vMoment");
    AuxillaryVariables["V"].push_back("HMoment");

    AuxillaryVariables["C"].resize(0);
    AuxillaryVariables["C"].push_back("Char1");
    AuxillaryVariables["C"].push_back("Char2");
    AuxillaryVariables["C"].push_back("Char3");
    AuxillaryVariables["C"].push_back("Char4");

    AuxillaryVariables["M"].resize(0);
    AuxillaryVariables["M"].push_back("uMoment");
    AuxillaryVariables["M"].push_back("vMoment");
    AuxillaryVariables["M"].push_back("qMoment");
    AuxillaryVariables["M"].push_back("PMoment");
    AuxillaryVariables["M"].push_back("HMoment");

  }

  return ;
}

void EulerSolver::RunLimiter() {
  if ( Limiter == "LiliaMoment") {
    field->resetPositivity(true);
    
    Run_LiliaMomentLimiter(DVx);
    Run_LiliaMomentLimiter(DVy);
    Run_LiliaMomentLimiter(DE);
    Run_LiliaMomentLimiter(D);
    //Run_LiliaMomentLimiter(T); // If needed, else compute it later using q and P ..
  }

  return ;
}

void EulerSolver::Run_LiliaMomentLimiter(int v) {
  field->computeMoments(v, Moment, CellMarker);
  //field->computeMoments(v, ModMoment);
  field->setFunctionsForVariables(1.0, Moment, Copy, ModMoment);
  field->limitMoments(Moment, ModMoment, CellMarker, (N+1)*(N+1)-1);
  field->convertMomentToVariable(ModMoment, v, CellMarker);

  return ;
}

void EulerSolver::RunPositivityLimiter() {
  
  if ( Limiter == "LiliaMoment") {
    field->ResetVariables_CellCentered(CellMarker, 1.5);

    updatePrimitiveVariables();
    field->resetPositivity(false);
    checkPositivity();
    Run_PositivityMomentLimiter(DVx, N+2);
    Run_PositivityMomentLimiter(DVy, N+2);
    Run_PositivityMomentLimiter(DE, N+2);
    Run_PositivityMomentLimiter(D, N+2);
    

    updatePrimitiveVariables();
    field->resetPositivity(false);
    checkPositivity();
    Run_PositivityMomentLimiter(DVx, 0);
    Run_PositivityMomentLimiter(DVy, 0);
    Run_PositivityMomentLimiter(DE, 0);
    Run_PositivityMomentLimiter(D, 0);
    //Run_PositivityMomentLimiter(T, N+2); // If needed, else compute it later using q and P ..
  }

  else if(Limiter == "CharacteristicLimiter") {
  
    field->computeMoments(DVx, uMoment, CellMarker);
    field->computeMoments(DVy, vMoment, CellMarker);
    field->computeMoments(D, qMoment, CellMarker);
    field->computeMoments(P, PMoment, CellMarker);
    field->computeMoments(DE, HMoment, CellMarker);

    // Finding gradient of  Density //Pressure
    field->delByDelX(P, dPdx, P, "central", Vx_plus_C);
    field->delByDelY(P, dPdy, P, "central", Vy_plus_C);
    field->computeMoments(dPdx, dPdxMoment, CellMarker);
    field->computeMoments(dPdy, dPdyMoment, CellMarker);


    field->findEigenMatrices(AuxillaryVariables["M"], CellMarker);

    // Find Characteristic Variables 
    field->convertVariabletoCharacteristic(AuxillaryVariables["V"], "Char1", 0, CellMarker);
    field->convertVariabletoCharacteristic(AuxillaryVariables["V"], "Char2", 4, CellMarker);
    field->convertVariabletoCharacteristic(AuxillaryVariables["V"], "Char3", 8, CellMarker);
    field->convertVariabletoCharacteristic(AuxillaryVariables["V"], "Char4", 12, CellMarker);

    // Limiting Characteristic Variables ;
    field->limitMoments(AuxillaryVariables["V"], Char1, CellMarker, 0);
    field->limitMoments(AuxillaryVariables["V"], Char2, CellMarker, 1);
    field->limitMoments(AuxillaryVariables["V"], Char3, CellMarker, 2);
    field->limitMoments(AuxillaryVariables["V"], Char4, CellMarker, 3);

    // update Conservative Variables
    field->convertCharacteristictoVariable(AuxillaryVariables["C"], "qMoment", 0, CellMarker);
    field->convertCharacteristictoVariable(AuxillaryVariables["C"], "uMoment", 4, CellMarker);
    field->convertCharacteristictoVariable(AuxillaryVariables["C"], "vMoment", 8, CellMarker);
    field->convertCharacteristictoVariable(AuxillaryVariables["C"], "HMoment", 12, CellMarker);

    field->convertMomentToVariable(qMoment, D, CellMarker);
    field->convertMomentToVariable(uMoment, DVx, CellMarker);
    field->convertMomentToVariable(vMoment, DVy, CellMarker);
    field->convertMomentToVariable(HMoment, DE, CellMarker);
  }


  return ;
}

void EulerSolver::Run_PositivityMomentLimiter(int v, unsigned Index) {
  field->computeMoments(v, Moment, CellMarker);
  field->scal(0.0, ModMoment);
  field->limitMoments(Moment, ModMoment, CellMarker, Index);
  field->convertMomentToVariable(ModMoment, v, CellMarker);

  return ;
}


void EulerSolver::FindL2Norm(function<double(double, double)> Density, function<double(double, double)> U ) {
  field->scal(0.0,ZERO);

  field->initializeVariable(DAnalytical, Density);
  field->initializeVariable(VxAnalytical, U);

  // Computing L2Norm

  double qNorm = 0.0, qAnalytic = 0.0;
  double uNorm = 0.0, uAnalytic = 0.0;
  
  qNorm = field->l2Norm(D, DAnalytical);
  qAnalytic = field->l2Norm(DAnalytical, ZERO);
  uNorm = field->l2Norm(Vx, VxAnalytical);
  uAnalytic = field->l2Norm(VxAnalytical, ZERO);

  cout << "L2norms :\n  Density : " <<qNorm <<" : L2 Norm : "<< qNorm/qAnalytic <<"\n";
  cout <<" u Velocity : " <<uNorm << " : L2 Norm :  " << uNorm/uAnalytic <<"\n";
  return ;
}

void EulerSolver::checkPositivity() {
  field->checkPositivity(D, CellMarker, "One");
  field->checkPositivity(P, CellMarker, "Two");

}