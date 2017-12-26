#include "../../includes/Solvers/EulerSolver.h"
#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/MaterialProperties.h"
#include "../../includes/Utilities/MathOperators.h"
#include "../../includes/Utilities/Thermodynamics.h"

// Primitve Variables
#define D   0
#define Vx  1
#define Vy  2
#define P   3
#define T   4

// Conservative Variables
#define DVx 5
#define DVy 6
#define DE  7
#define De  8
#define KE  9

// Inviscid Flux
#define DVxVx_plus_P 10
#define DVxVy        11
#define DVyVy_plus_P 12
#define DE_plus_P_Vx 13
#define DE_plus_P_Vy 14

// Eigen Values
#define C         15
#define Vx_plus_C 16
#define Vy_plus_C 17

// Aux. Variables // Define all auxillary (needed or not here)
#define K1D   18
#define K1DVx 19
#define K1DVy 20
#define K1DE  21
#define K2D   22
#define K2DVx 23
#define K2DVy 24
#define K2DE  25
#define K3D   26
#define K3DVx 27
#define K3DVy 28
#define K3DE  29
#define dbydx 30
#define dbydy 31

#define DAnalytical  32
#define VxAnalytical 33
#define ZERO         34

// Shock Detector
#define CellMarkerG  35

// Limiter
#define Moment    36
#define ModMoment 37
// Add corresponding variables for Characteristic Limiter

// cell centered variables 
#define Dx         0 
#define Dt         1
#define UMax       2
#define CellMarker 3




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
  field->addVariable_withBounary(D);
  field->addVariable_withBounary(Vx);
  field->addVariable_withBounary(Vy);
  field->addVariable_withBounary(P);
  field->addVariable_withBounary(T);
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
  field->addVariable_withBounary(DVx);
  field->addVariable_withBounary(DVy);
  field->addVariable_withBounary(DE);
  field->addVariable_withBounary(De); // Internal Energy, added just for ease of manipulating Energy, Remove later if not required !!
  field->addVariable_withBounary(KE); // Kinetic Energy , " "

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
  field->addVariable_CellCentered(Dx);
  field->FindMindx(Dx);
  field->addVariable_CellCentered(Dt);
  field->addVariable_CellCentered(UMax);
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
  field->setFunctionsForVariables(1.0, KE, 1.0, De, Add, DE);
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
  setInternalEnergy();
  setEnergy();
  return ;
}

void EulerSolver::setInviscidFlux() {
  //cout << " Calling setInviscidFlux " << endl;
  field->addVariable_withBounary(DVxVx_plus_P);
  field->addVariable_withBounary(DVxVy);
  field->addVariable_withBounary(DVyVy_plus_P);
  field->addVariable_withBounary(DE_plus_P_Vx);
  field->addVariable_withBounary(DE_plus_P_Vy);

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
  field->addVariable_withoutBounary(K1D);
  field->addVariable_withoutBounary(K1DVx);
  field->addVariable_withoutBounary(K1DVy);
  field->addVariable_withoutBounary(K1DE);
  field->addVariable_withoutBounary(K2D);
  field->addVariable_withoutBounary(K2DVx);
  field->addVariable_withoutBounary(K2DVy);
  field->addVariable_withoutBounary(K2DE);
  field->addVariable_withoutBounary(K3D);
  field->addVariable_withoutBounary(K3DVx);
  field->addVariable_withoutBounary(K3DVy);
  field->addVariable_withoutBounary(K3DE);
  field->addVariable_withoutBounary(dbydx);
  field->addVariable_withoutBounary(dbydy);

  field->addVariable_withBounary(DAnalytical);
  field->addVariable_withBounary(VxAnalytical);
  
  field->addVariable_withoutBounary(ZERO);

  //field->addVariable_withBounary("xFlux");
  //field->addVariable_withBounary("yFlux");
  return ;
}

void EulerSolver::setEigenValues() {
 // cout << "Calling setEigenValues " << endl;
  field->addVariable_withBounary(C);
  field->addVariable_withBounary(Vx_plus_C);
  field->addVariable_withBounary(Vy_plus_C);// Recheck formulation of eigen value !!

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

  //field->scal(0.0, K);
  //field->axpy(-1.0, dbydx, K);
  //field->axpy(-1.0, dbydy, K);
  field->setFunctionsForVariables(-1.0, dbydx, -1.0, dbydy, Add, K);
  
  //field->axpy(0.5*dt, K, Var);
  field->setFunctionsForVariables(0.5*dt, K, 1.0, Var, Add, Var);

  return;
}

void EulerSolver::RK_Step2(int Var, int FluxX, int FluxY, int K1, int K2) {
  field->delByDelX(FluxX, dbydx, Var, "rusanov", Vx_plus_C);
  field->delByDelY(FluxY, dbydy, Var, "rusanov", Vy_plus_C);

  //field->scal(0.0, K2);
  //field->axpy(-1.0, dbydx, K2);
  //field->axpy(-1.0, dbydy, K2);
  field->setFunctionsForVariables(-1.0, dbydx, -1.0, dbydy, Add, K2);
  
  //field->axpy(-1.5*dt, K1, Var);
  //field->axpy(2.0*dt, K2, Var);
  field->setFunctionsForVariables(-1.5*dt, K1, 2.0*dt, K2, 1.0, Var, Add, Var);

  return;
}

void EulerSolver::RK_Step3(int Var, int FluxX, int FluxY, int K1, int K2, int K3) {
  field->delByDelX(FluxX, dbydx, Var, "rusanov", Vx_plus_C);
  field->delByDelY(FluxY, dbydy, Var, "rusanov", Vy_plus_C);

  //field->scal(0.0, K3);
  //field->axpy(-1.0, dbydx, K3);
  //field->axpy(-1.0, dbydy, K3);
  field->setFunctionsForVariables(-1.0, dbydx, -1.0, dbydy, Add, K3);
  
  //field->axpy((7.0/6.0)*dt, K1, Var);
  //field->axpy(-(4.0/3.0)*dt, K2, Var);
  //field->axpy((1.0/6.0)*dt, K3, Var);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1, -(4.0/3.0)*dt, K2, (1.0/6.0)*dt, K3, 1.0, Var, Add, Var);
  return;
}

void EulerSolver::setTimeStep() {
  field->setFunctionsForCellCenterVariablesfromDomainVariables(1.0, Vx_plus_C, 1.0, Vy_plus_C, Maximum, UMax);
  field->FindTimestep(Dt, Dx, UMax, CFL);
  dt = field->FindMindt(Dt);
  //dt = CFL * dx/speed;
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
    field->addVariable_CellCentered(CellMarker);
    field->addVariable_withBounary(CellMarkerG);
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

  field->updateOutFlowBoundary(Vx, Vy);
  field->updateCellMarker(D, CellMarker);

  return ;
}



void EulerSolver::SetLimiter(string _Limiter) {
  Limiter = _Limiter;
  SetLimiterVariables();
  return ;
}

void EulerSolver::SetLimiterVariables() {
  if (Limiter == "LiliaMoment") {
    /*field->addVariable_CellCentered(CellMarker);
    field->ResetVariables_CellCentered(CellMarker, 1.5);
    field->addVariable_withBounary(CellMarkerG);
    field->scal(0.0, CellMarkerG);*/
    //field->addVariable_CellCentered("Max");
    //field->addVariable_CellCentered("Min");


    field->addVariable_withBounary(Moment);
    field->addVariable_withBounary(ModMoment);
    field->setVanderMandMatrix();
  }

  return ;
}

void EulerSolver::RunLimiter() {
  if ( Limiter == "LiliaMoment") {
    field->resetPositivity();
    Run_LiliaMomentLimiter(DVx);
    Run_LiliaMomentLimiter(DVy);
    Run_LiliaMomentLimiter(DE);
    Run_LiliaMomentLimiter(D);
    //Run_LiliaMomentLimiter(T); // If needed, else compute it later using q and P ..
  }

  return ;
}

void EulerSolver::Run_LiliaMomentLimiter(int v) {
  field->computeMoments(v, Moment);
  field->computeMoments(v, ModMoment);
  field->limitMoments(Moment, ModMoment, CellMarker, (N+1)*(N+1)-1);
  field->convertMomentToVariable(ModMoment, v, CellMarker);

  return ;
}

void EulerSolver::RunPositivityLimiter() {
  field->resetPositivity();
  if ( Limiter == "LiliaMoment") {
    updatePrimitiveVariables();
    checkPositivity();
    Run_PositivityMomentLimiter(DVx, N+2);
    Run_PositivityMomentLimiter(DVy, N+2);
    Run_PositivityMomentLimiter(DE, N+2);
    Run_PositivityMomentLimiter(D, N+2);

    updatePrimitiveVariables();
    checkPositivity();
    Run_PositivityMomentLimiter(DVx, 0);
    Run_PositivityMomentLimiter(DVy, 0);
    Run_PositivityMomentLimiter(DE, 0);
    Run_PositivityMomentLimiter(D, 0);
    //Run_PositivityMomentLimiter(T, N+2); // If needed, else compute it later using q and P ..
  }

  return ;
}

void EulerSolver::Run_PositivityMomentLimiter(int v, unsigned Index) {
  field->computeMoments(v, Moment);
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