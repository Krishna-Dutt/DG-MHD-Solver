#include "../../includes/Solvers/NavierStokesSolver.h"
#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/MaterialProperties.h"
#include "../../includes/Utilities/MathOperators.h"
#include "../../includes/Utilities/ThermodynamicFunctions.h"


NSSolver::NSSolver(int _ne_x, int _ne_y, int _N) {
    ne_x = _ne_x;
    ne_y = _ne_y;
    N = _N;
    time = 0.0;
}

NSSolver::~NSSolver() {
  delete field;
}

void NSSolver::setDomain(double _x1, double _y1, double _x2, double _y2) {
    x1 = _x1;
    y1 = _y1;
    x2 = _x2;
    y2 = _y2;
    field = new DG_Field_2d(ne_x, ne_y, N, x1, y1, x2, y2);
    field->setSystem("EULER");
    Dimension = 4;

   return ;
}

void NSSolver::setPrimitiveVariables(){
  D  = field->addVariable_withBounary("q");
  Vx = field->addVariable_withBounary("u");
  Vy = field->addVariable_withBounary("v");
  P  = field->addVariable_withBounary("P");
  T  = field->addVariable_withBounary("T");
  return ;
}

void NSSolver::setGradientPrimitiveVariables(){
  dDdx = field->addVariable_withBounary("dqdx");
  dVxdx = field->addVariable_withBounary("dudx");
  dVydx = field->addVariable_withBounary("dvdx");
  dPdx = field->addVariable_withBounary("dPdx");
  //field->addVariable_withoutBounary("dTdx");

  dDdy = field->addVariable_withBounary("dqdy");
  dVxdy = field->addVariable_withBounary("dudy");
  dVydy = field->addVariable_withBounary("dvdy");
  dPdy = field->addVariable_withBounary("dPdy");
  //field->addVariable_withoutBounary("dTdy");

  return ;
}


void NSSolver::setConservativeVariables(){
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

/*void NSSolver::setMaterialPropertyVariables(){
  field->addVariable_withoutBounary("meu");
  return ;
}*/

void NSSolver::setInitialDensity(function<double(double, double)> Rho) {
    field->initializeVariable(D, Rho);
    return ;
}

void NSSolver::setInitialPressure(function<double(double, double)> Pr) {
    field->initializeVariable(P, Pr);
    return ;
}

void NSSolver::setInitialTemperature(function<double(double, double)> Tp) {
    field->initializeVariable(T, Tp);
    return ;
}

void NSSolver::setInitialVelocity(function<double(double, double)> U, function<double(double, double)> V) {
    field->initializeVariable(Vx, U);
    field->initializeVariable(Vy, V);
    return ;
}

void NSSolver::setBoundaryCondtions(string type1, string type2, string type3, string type4) {
      // Setting BC as Outflow type to test Methods
    field->setBottomBoundary(type1);
    field->setRightBoundary(type2);
    field->setTopBoundary(type3);
    field->setLeftBoundary(type4);

    return ;
}

void NSSolver::setSolver(double _CFL, double _time, int _no_of_time_steps) {
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

void NSSolver::setXMomentum() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vx, Product, DVx);
  return ;
}

void NSSolver::setYMomentum() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vy, Product, DVy);
  return ;
}

void NSSolver::setEnergy() {
  field->setFunctionsForVariables(1.0, KE, 1.0, De, Addab, DE);
  return ;
}

void NSSolver::setInternalEnergyfromPrimitive() {
  field->setFunctionsForVariables(1.0, D, 1.0, T, 1.0, P, IE, De);
  return ;
}

void NSSolver::setInternalEnergy() {
  field->setFunctionsForVariables(1.0, DE, 1.0, KE, Subtract, De);
  return ;
}


void NSSolver::setKineticEnergy() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vx, 1.0, Vy, KineticEnergy, KE);
  return ;
}

void NSSolver::updateVelocity() {
  field->setFunctionsForVariables(1.0, DVx, 1.0, D, Divide, Vx);
  field->setFunctionsForVariables(1.0 ,DVy, 1.0, D, Divide, Vy);
  return ;
}

void NSSolver::updateTemperature() {
  field->setFunctionsForVariables(1.0, De, 1.0, D, Temperature, T);
  return ;
}

void NSSolver::updatePressure() {
  field->setFunctionsForVariables(1.0, D, 1.0, De, Pressure, P);
  return ;
}

void NSSolver::updatePrimitiveVariables() {
  updateVelocity();
  setKineticEnergy();
  setInternalEnergy();
  updateTemperature();
  updatePressure();
  return ;
}

void NSSolver::updateConservativeVariables() {
  setXMomentum();
  setYMomentum();
  setKineticEnergy();
  setInternalEnergyfromPrimitive();
  setEnergy();
  return ;
}

void NSSolver::setInviscidFlux() {
  //cout << " Calling setInviscidFlux " << endl;
  DVxVx_plus_P = field->addVariable_withBounary("quu_plus_P");
  DVxVy        = field->addVariable_withBounary("quv");
  DVyVy_plus_P = field->addVariable_withBounary("qvv_plus_P");
  DE_plus_P_Vx = field->addVariable_withBounary("qE_plus_P_u");
  DE_plus_P_Vy = field->addVariable_withBounary("qE_plus_P_v");

  //updateInviscidFlux();
  return ;
}

void NSSolver::updateInviscidFlux() {
  field->setFunctionsForVariables(1.0, DVx, 1.0, Vx, 1.0, P, MomentumFluxPressure, DVxVx_plus_P);
  field->setFunctionsForVariables(1.0, DVy, 1.0, Vy, 1.0, P, MomentumFluxPressure, DVyVy_plus_P);
  field->setFunctionsForVariables(1.0, DVx, 1.0, Vy, MomentumFlux, DVxVy);
  field->setFunctionsForVariables(1.0, DE, 1.0, P, 1.0, Vx, EnergyFlux, DE_plus_P_Vx);
  field->setFunctionsForVariables(1.0, DE, 1.0, P, 1.0, Vy, EnergyFlux, DE_plus_P_Vy);
  
  return ;
}

void NSSolver::setViscousFlux() {
  TauXX = field->addVariable_withBounary("Tauxx");
  TauXY = field->addVariable_withBounary("Tauxy");
  TauYY = field->addVariable_withBounary("Tauyy");
  EViscX = field->addVariable_withBounary("Eviscousx");
  EViscY = field->addVariable_withBounary("Eviscousy");

  //updateInviscidFlux();
  return ;
}

void NSSolver::updateViscousFlux() {
  field->setFunctionsForVariables((4.0/3.0), dVxdx, (-2.0/3.0), dVydy, ViscousStress, TauXX);
  field->setFunctionsForVariables((4.0/3.0), dVydy, (-2.0/3.0), dVxdx, ViscousStress, TauYY);
  field->setFunctionsForVariables(1.0, dVxdy, 1.0, dVydx, ViscousStress, TauXY);
  field->setFunctionsForVariables(1.0, Vx, 1.0, TauXX, 1.0, Vy, 1.0, TauXY, EnergyViscous, EViscX);
  field->setFunctionsForVariables(1.0, Vx, 1.0, TauXY, 1.0, Vy, 1.0, TauYY, EnergyViscous, EViscY);
  
  return ;
}

void NSSolver::updatePrimitiveGradient() {
  field->delByDelX(D, dDdx, "central");
  field->delByDelX(Vx, dVxdx, "central");
  field->delByDelX(Vy, dVydx, "central");
  field->delByDelX(P, dPdx, "central");
  //field->delByDelX(T, "dTdx", "central");

  field->delByDelY(D, dDdy, "central");
  field->delByDelY(Vx, dVxdy, "central");
  field->delByDelY(Vy, dVydy, "central");
  field->delByDelY(P, dPdy, "central");
  //field->delByDelY(T, "dTdy", "central");

  return ;
}

void NSSolver::setAuxillaryVariables() {
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
  dbydxD = field->addVariable_withoutBounary();
  dbydyD = field->addVariable_withoutBounary();
  dbydxDVx = field->addVariable_withoutBounary();
  dbydyDVx = field->addVariable_withoutBounary();
  dbydxDVy = field->addVariable_withoutBounary();
  dbydyDVy = field->addVariable_withoutBounary();
  dbydxDE = field->addVariable_withoutBounary();
  dbydyDE = field->addVariable_withoutBounary();

  DAnalytical  = field->addVariable_withBounary("qAnalytic");
  VxAnalytical = field->addVariable_withBounary("uAnalytic");
  
  ZERO = field->addVariable_withoutBounary();
  
  return ;
}

void NSSolver::setEigenValues() {
 // cout << "Calling setEigenValues " << endl;
  C         = field->addVariable_withBounary("c");
  Vx_plus_C = field->addVariable_withBounary("u_plus_c");
  Vy_plus_C = field->addVariable_withBounary("v_plus_c");// Recheck formulation of eigen value !!

  updateEigenValues();
  return ;
}

void NSSolver::updateEigenValues() {
 // cout << " Calling updateEigenValues " << endl;
  field->setFunctionsForVariables(1.0, D, 1.0, P, SoundSpeed, C);
  field->setFunctionsForVariables(1.0, Vx, 1.0, C, ModulusAdd, Vx_plus_C);
  field->setFunctionsForVariables(1.0, Vy, 1.0, C, ModulusAdd, Vy_plus_C);
  return ;
}

void NSSolver::RK_Step1() {

  // Inviscid Flux
  int FluxX[] = { DVxVx_plus_P, DVxVy, DE_plus_P_Vx, DVx};
  int FluxY[] = { DVxVy, DVyVy_plus_P, DE_plus_P_Vy, DVy};
  int FluxVarx[] = {Vx_plus_C};
  int FluxVary[] = {Vy_plus_C};
  int Var[] = { DVx, DVy, DE, D};
  int DbyDx[] = { dbydxDVx, dbydxDVy, dbydxDE, dbydxD};
  int DbyDy[] = { dbydyDVx, dbydyDVy, dbydyDE, dbydyD};

  field->delByDelX(FluxX, DbyDx, Var, "rusanov", FluxVarx, 4);
  field->delByDelY(FluxY, DbyDy, Var, "rusanov", FluxVary, 4);

  field->setFunctionsForVariables(-1.0, dbydxD, -1.0, dbydyD, Addab, K1D);
  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, Addab, K1DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, Addab, K1DVy);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, Addab, K1DE);

  // Viscous Flux
  int VFluxX[] = { TauXX, TauXY, EViscX};
  int VFluxY[] = { TauXY, TauYY, EViscY};

  field->delByDelX(VFluxX, DbyDx, Var, "central", FluxVarx, 3);
  field->delByDelY(VFluxY, DbyDy, Var, "central", FluxVary, 3);

  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, 1.0, K1DVx, Addabc, K1DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, 1.0, K1DVy, Addabc, K1DVy);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, 1.0, K1DE, Addabc, K1DE);


  
  field->setFunctionsForVariables(0.5*dt, K1D, 1.0, D, Addab, D);
  field->setFunctionsForVariables(0.5*dt, K1DVx, 1.0, DVx, Addab, DVx);
  field->setFunctionsForVariables(0.5*dt, K1DVy, 1.0, DVy, Addab, DVy);
  field->setFunctionsForVariables(0.5*dt, K1DE, 1.0, DE, Addab, DE);

  return;
}

void NSSolver::RK_Step2() {
  // Inviscid Flux
  int FluxX[] = { DVxVx_plus_P, DVxVy, DE_plus_P_Vx, DVx};
  int FluxY[] = { DVxVy, DVyVy_plus_P, DE_plus_P_Vy, DVy};
  int FluxVarx[] = {Vx_plus_C};
  int FluxVary[] = {Vy_plus_C};
  int Var[] = { DVx, DVy, DE, D};
  int DbyDx[] = { dbydxDVx, dbydxDVy, dbydxDE, dbydxD};
  int DbyDy[] = { dbydyDVx, dbydyDVy, dbydyDE, dbydyD};

  field->delByDelX(FluxX, DbyDx, Var, "rusanov", FluxVarx, 4);
  field->delByDelY(FluxY, DbyDy, Var, "rusanov", FluxVary, 4);
  
  field->setFunctionsForVariables(-1.0, dbydxD, -1.0, dbydyD, Addab, K2D);
  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, Addab, K2DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, Addab, K2DVy);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, Addab, K2DE);

  // Viscous Flux
  int VFluxX[] = { TauXX, TauXY, EViscX};
  int VFluxY[] = { TauXY, TauYY, EViscY};

  field->delByDelX(VFluxX, DbyDx, Var, "central", FluxVarx, 3);
  field->delByDelY(VFluxY, DbyDy, Var, "central", FluxVary, 3);

  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, 1.0, K2DVx, Addabc, K2DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, 1.0, K2DVy, Addabc, K2DVy);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, 1.0, K2DE, Addabc, K2DE);

  
  field->setFunctionsForVariables(-1.5*dt, K1D, 2.0*dt, K2D, 1.0, D, Addabc, D);
  field->setFunctionsForVariables(-1.5*dt, K1DVx, 2.0*dt, K2DVx, 1.0, DVx, Addabc, DVx);
  field->setFunctionsForVariables(-1.5*dt, K1DVy, 2.0*dt, K2DVy, 1.0, DVy, Addabc, DVy);
  field->setFunctionsForVariables(-1.5*dt, K1DE, 2.0*dt, K2DE, 1.0, DE, Addabc, DE);

  return;
}

void NSSolver::RK_Step3() {
  // Inviscid Flux
  int FluxX[] = { DVxVx_plus_P, DVxVy, DE_plus_P_Vx, DVx};
  int FluxY[] = { DVxVy, DVyVy_plus_P, DE_plus_P_Vy, DVy};
  int FluxVarx[] = {Vx_plus_C};
  int FluxVary[] = {Vy_plus_C};
  int Var[] = { DVx, DVy, DE, D};
  int DbyDx[] = { dbydxDVx, dbydxDVy, dbydxDE, dbydxD};
  int DbyDy[] = { dbydyDVx, dbydyDVy, dbydyDE, dbydyD};

  field->delByDelX(FluxX, DbyDx, Var, "rusanov", FluxVarx, 4);
  field->delByDelY(FluxY, DbyDy, Var, "rusanov", FluxVary, 4);

  field->setFunctionsForVariables(-1.0, dbydxD, -1.0, dbydyD, Addab, K3D);
  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, Addab, K3DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, Addab, K3DVy);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, Addab, K3DE);

  // Viscous Flux
  int VFluxX[] = { TauXX, TauXY, EViscX};
  int VFluxY[] = { TauXY, TauYY, EViscY};

  field->delByDelX(VFluxX, DbyDx, Var, "central", FluxVarx, 3);
  field->delByDelY(VFluxY, DbyDy, Var, "central", FluxVary, 3);

  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, 1.0, K3DVx, Addabc, K3DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, 1.0, K3DVy, Addabc, K3DVy);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, 1.0, K3DE, Addabc, K3DE);

  
  field->setFunctionsForVariables((7.0/6.0)*dt, K1D, -(4.0/3.0)*dt, K2D, (1.0/6.0)*dt, K3D, 1.0, D, Addabcd, D);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1DVx, -(4.0/3.0)*dt, K2DVx, (1.0/6.0)*dt, K3DVx, 1.0, DVx, Addabcd, DVx);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1DVy, -(4.0/3.0)*dt, K2DVy, (1.0/6.0)*dt, K3DVy, 1.0, DVy, Addabcd, DVy);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1DE, -(4.0/3.0)*dt, K2DE, (1.0/6.0)*dt, K3DE, 1.0, DE, Addabcd, DE);

  return;
}

void NSSolver::setTimeStep() {
  field->setFunctionsForCellCenterVariablesfromDomainVariables(1.0, Vx_plus_C, 1.0, Vy_plus_C, Maximum, UMax);
  field->FindTimestep(Dt, Dx, UMax, CFL);
  dt = field->FindMindt(Dt);

  return ;
}


void NSSolver::solve() {
  // Requires all Primitive and Conservative Variables to be setup and initialised.
  double t = 0.0;
  int count = 0;
  
 
  // Till Now all variables have to be initialised !!
  // For loop to march in time !!
  while(t <= time) {
    // First Step of RK3

    updateInviscidFlux();
    updatePrimitiveGradient(); 
    updateViscousFlux();
    updateEigenValues();

    if ( count%no_of_time_steps == 0) {
      setTimeStep();
      cout << "Time Step : " << dt << " , Time : " << t << "\n"; 
      count = 0;
    }

    RK_Step1();
    
    RunShockDetector();
    RunLimiter();
    RunPositivityLimiter();

    field->updateBoundary(t);
    updatePrimitiveVariables();
   
    
    // Second Step of RK3
    updateInviscidFlux();
    updatePrimitiveGradient(); 
    updateViscousFlux();
    updateEigenValues();
    
    RK_Step2();
    
    RunShockDetector();
    RunLimiter();
    RunPositivityLimiter(); 
    
    field->updateBoundary(t);
    updatePrimitiveVariables();

   // Third (Final) Step of RK3
    updateInviscidFlux();
    updatePrimitiveGradient(); 
    updateViscousFlux();
    updateEigenValues();

    RK_Step3();
    
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

void NSSolver::plot(string filename) {
    field->writeVTK(filename);
    return ;
}

void NSSolver::SetShockDetector(string _ShockDetector) {
  ShockDetector = _ShockDetector;
  SetShockDetectorVariables();
  return ;
}

void NSSolver::SetShockDetectorVariables() {
  
  if (ShockDetector == "KXRCF") {
    CellMarker  = field->addVariable_CellCentered();
    CellMarkerG = field->addVariable_withBounary("cellmarker");

    field->scal(0.0, CellMarkerG);
    
  }

  return ;
}

void NSSolver::RunShockDetector() {
  if (ShockDetector == "KXRCF") {
    Run_KXRCF();
 }

  return ;
}

void NSSolver::Run_KXRCF() {
 
  field->ResetVariables_CellCentered(CellMarker, 1.5);
  field->ResetMap_OutFlow();

  //field->updateOutFlowBoundary(Vx, Vy);
  //field->updateCellMarker(D, CellMarker);

  return ;
}



void NSSolver::SetLimiter(string _Limiter) {
  Limiter = _Limiter;
  SetLimiterVariables();
  return ;
}

void NSSolver::SetLimiterVariables() {
  if (Limiter == "LiliaMoment") {
    /*CellMarker  = field->addVariable_CellCentered();
    field->ResetVariables_CellCentered(CellMarker, 1.5);
    CellMarkerG = field->addVariable_withBounary("cellmarker");
    field->scal(0.0, CellMarkerG);*/

    Moment    = field->addVariable_withBounary("moment");
    field->scal(0.0, Moment);
    ModMoment = field->addVariable_withBounary("modmoment");
    field->scal(0.0, ModMoment);
    uMoment = field->addVariable_withoutBounary();
    vMoment = field->addVariable_withoutBounary();
    qMoment = field->addVariable_withoutBounary();
    HMoment = field->addVariable_withoutBounary();
    uModMoment = field->addVariable_withoutBounary();
    vModMoment = field->addVariable_withoutBounary();
    qModMoment = field->addVariable_withoutBounary();
    HModMoment = field->addVariable_withoutBounary();
    field->scal(0.0, uMoment);
    field->scal(0.0, vMoment);
    field->scal(0.0, qMoment);
    field->scal(0.0, HMoment);
    field->scal(0.0, uModMoment);
    field->scal(0.0, vModMoment);
    field->scal(0.0, qModMoment);
    field->scal(0.0, HModMoment);
    field->setVanderMandMatrix();
  }
  else if (Limiter == "CharacteristicLimiter") {
    field->setEigenMatrices(Dimension);
    field->setVanderMandMatrix();
    Moment    = field->addVariable_withBounary("moment");
    field->scal(0.0, Moment);
    ModMoment = field->addVariable_withBounary("modmoment");
    field->scal(0.0, ModMoment);
    // Change later by adding separate settings for positivity limiter!!
    uMoment = field->addVariable_withoutBounary();
    vMoment = field->addVariable_withoutBounary();
    qMoment = field->addVariable_withoutBounary();
    HMoment = field->addVariable_withoutBounary();
    field->scal(0.0, uMoment);
    field->scal(0.0, vMoment);
    field->scal(0.0, qMoment);
    field->scal(0.0, HMoment);
    uModMoment = field->addVariable_withoutBounary();
    vModMoment = field->addVariable_withoutBounary();
    qModMoment = field->addVariable_withoutBounary();
    HModMoment = field->addVariable_withoutBounary();
    field->scal(0.0, uModMoment);
    field->scal(0.0, vModMoment);
    field->scal(0.0, qModMoment);
    field->scal(0.0, HModMoment);
    

    Char1 = field->addVariable_withoutBounary();
    Char2 = field->addVariable_withoutBounary();
    Char3 = field->addVariable_withoutBounary();
    Char4 = field->addVariable_withoutBounary();
    field->scal(0.0, Char1);
    field->scal(0.0, Char2);
    field->scal(0.0, Char3);
    field->scal(0.0, Char4);

    //dPdx = field->addVariable_withBounary("dPdx");
    //dPdy = field->addVariable_withBounary("dPdy");
    dPdxMoment = field->addVariable_withoutBounary();
    dPdyMoment = field->addVariable_withoutBounary();
    field->scal(0.0, dPdxMoment);
    field->scal(0.0, dPdyMoment);

  }

  return ;
}

void NSSolver::RunLimiter() {
  field->resetPositivity(true);
  
  if ( Limiter == "LiliaMoment") {
        
    /*Run_LiliaMomentLimiter(DVx);
    Run_LiliaMomentLimiter(DVy);
    Run_LiliaMomentLimiter(DE);
    Run_LiliaMomentLimiter(D);
    */
    //Run_LiliaMomentLimiter(T); // If needed, else compute it later using q and P ..
    int Var[] ={D, DVx, DVy, DE};
    int Mom[] = {qMoment, uMoment, vMoment, HMoment};
    int ModMom[] = {qModMoment, uModMoment, vModMoment, HModMoment}; 
    field->computeMoments(Var, Mom, CellMarker, 4);
    //field->computeMoments(v, ModMoment);
    field->setFunctionsForVariables(1.0, qMoment, Copy, qModMoment);
    field->setFunctionsForVariables(1.0, uMoment, Copy, uModMoment);
    field->setFunctionsForVariables(1.0, vMoment, Copy, vModMoment);
    field->setFunctionsForVariables(1.0, HMoment, Copy, HModMoment);

    field->limitMoments(Mom, ModMom, CellMarker, (N+1)*(N+1)-1, 4);
    field->convertMomentToVariable(ModMom, Var, CellMarker, 4);

  }

  else if(Limiter == "CharacteristicLimiter") {
    int AuxVarM[6] = {uMoment, vMoment, qMoment, HMoment, dPdxMoment, dPdyMoment};
    int AuxVarV[4] = {qMoment, uMoment, vMoment, HMoment};
    int AuxVarC[4] = {Char1, Char2, Char3, Char4};
  
    field->computeMoments(DVx, uMoment, CellMarker);
    field->computeMoments(DVy, vMoment, CellMarker);
    field->computeMoments(D, qMoment, CellMarker);
    field->computeMoments(DE, HMoment, CellMarker);

    // Finding gradient of  Density //Pressure
    field->delByDelX(P, dPdx, P, "central", Vx_plus_C);
    field->delByDelY(P, dPdy, P, "central", Vy_plus_C);
    field->computeMoments(dPdx, dPdxMoment, CellMarker);
    field->computeMoments(dPdy, dPdyMoment, CellMarker);


    field->findEigenMatrices(AuxVarM, CellMarker);

    // Find Characteristic Variables 
    field->convertVariabletoCharacteristic(AuxVarV, AuxVarC, 0, CellMarker);
    
    // Limiting Characteristic Variables ;
    field->limitMoments(AuxVarV, AuxVarC, CellMarker, 0);
    
    // update Conservative Variables
    field->convertCharacteristictoVariable(AuxVarC, AuxVarV, 0, CellMarker);
    
    field->convertMomentToVariable(qMoment, D, CellMarker);
    field->convertMomentToVariable(uMoment, DVx, CellMarker);
    field->convertMomentToVariable(vMoment, DVy, CellMarker);
    field->convertMomentToVariable(HMoment, DE, CellMarker);
  }


  return ;
}

void NSSolver::Run_LiliaMomentLimiter(int v) {
  field->computeMoments(v, Moment, CellMarker);
  //field->computeMoments(v, ModMoment);
  field->setFunctionsForVariables(1.0, Moment, Copy, ModMoment);
  field->limitMoments(Moment, ModMoment, CellMarker, (N+1)*(N+1)-1);
  field->convertMomentToVariable(ModMoment, v, CellMarker);

  return ;
}

void NSSolver::RunPositivityLimiter() {
  
  if ( Limiter == "LiliaMoment" || Limiter == "CharacteristicLimiter") {
    int Var[] ={D, DVx, DVy, DE};
    int Mom[] = {qMoment, uMoment, vMoment, HMoment};
    int ModMom[] = {qModMoment, uModMoment, vModMoment, HModMoment}; 

    field->ResetVariables_CellCentered(CellMarker, 1.5);

    updatePrimitiveVariables();
    field->resetPositivity(false);
    checkPositivity();
    /*Run_PositivityMomentLimiter(DVx, N+2);
    Run_PositivityMomentLimiter(DVy, N+2);
    Run_PositivityMomentLimiter(DE, N+2);
    Run_PositivityMomentLimiter(D, N+2);*/

    field->computeMoments(Var, Mom, CellMarker, 4);
    field->scal(0.0, qModMoment);
    field->scal(0.0, uModMoment);
    field->scal(0.0, vModMoment);
    field->scal(0.0, HModMoment);
    field->limitMoments(Mom, ModMom, CellMarker, N+2, 4);
    field->convertMomentToVariable(ModMom, Var, CellMarker, 4);


    updatePrimitiveVariables();
    field->resetPositivity(false);
    checkPositivity();
    /*Run_PositivityMomentLimiter(DVx, 0);
    Run_PositivityMomentLimiter(DVy, 0);
    Run_PositivityMomentLimiter(DE, 0);
    Run_PositivityMomentLimiter(D, 0);*/
    //Run_PositivityMomentLimiter(T, N+2); // If needed, else compute it later using q and P ..
    field->computeMoments(Var, Mom, CellMarker, 4);
    field->scal(0.0, qModMoment);
    field->scal(0.0, uModMoment);
    field->scal(0.0, vModMoment);
    field->scal(0.0, HModMoment);
    field->limitMoments(Mom, ModMom, CellMarker, 0, 4);
    field->convertMomentToVariable(ModMom, Var, CellMarker, 4);
  }

  return ;
}

void NSSolver::Run_PositivityMomentLimiter(int v, unsigned Index) {
  field->computeMoments(v, Moment, CellMarker);
  field->scal(0.0, ModMoment);
  field->limitMoments(Moment, ModMoment, CellMarker, Index);
  field->convertMomentToVariable(ModMoment, v, CellMarker);

  return ;
}


void NSSolver::FindL2Norm(function<double(double, double)> Density, function<double(double, double)> U ) {
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

void NSSolver::checkPositivity() {
  field->checkPositivity(D, CellMarker, "One");
  field->checkPositivity(P, CellMarker, "Two");

}