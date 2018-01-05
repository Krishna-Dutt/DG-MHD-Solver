#include "../../includes/Solvers/IdealMHDSolver.h"
#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/MaterialProperties.h"
#include "../../includes/Utilities/MathOperators.h"
#include "../../includes/Utilities/ThermodynamicFunctions.h"


IdealMHDSolver::IdealMHDSolver(int _ne_x, int _ne_y, int _N) {
    ne_x = _ne_x;
    ne_y = _ne_y;
    N = _N;
    time = 0.0;
}

IdealMHDSolver::~IdealMHDSolver() {
  delete field;
}

void IdealMHDSolver::setDomain(double _x1, double _y1, double _x2, double _y2) {
    x1 = _x1;
    y1 = _y1;
    x2 = _x2;
    y2 = _y2;
    field = new DG_Field_2d(ne_x, ne_y, N, x1, y1, x2, y2);
    field->setSystem("EULER");//change later to MHD 
    Dimension = 4;

   return ;
}

void IdealMHDSolver::setPrimitiveVariables(){
  D  = field->addVariable_withBounary("q");
  Vx = field->addVariable_withBounary("u");
  Vy = field->addVariable_withBounary("v");
  P  = field->addVariable_withBounary("P");
  T  = field->addVariable_withBounary("T");
  Bx  = field->addVariable_withBounary("Bx");
  By  = field->addVariable_withBounary("By");
  Bz  = field->addVariable_withBounary("Bz");

  Pt  = field->addVariable_withBounary("Pt");
  
  return ;
}

/*void IdealMHDSolver::setGradientPrimitiveVariables(){
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


void IdealMHDSolver::setConservativeVariables(){
  DVx = field->addVariable_withBounary("qu");
  DVy = field->addVariable_withBounary("qv");
  DE  = field->addVariable_withBounary("qE");
  // Check definition of internal and total energy !!
  De  = field->addVariable_withBounary("qe"); // Internal Energy, added just for ease of manipulating Energy, Remove later if not required !!
  KE  = field->addVariable_withBounary("KE"); // Kinetic Energy , " "

  field->addConservativeVariables(DVx);
  field->addConservativeVariables(DVy);
  field->addConservativeVariables(DE);
  field->addConservativeVariables(Bx);
  field->addConservativeVariables(By);
  field->addConservativeVariables(Bz);

  return ;
}

/*void IdealMHDSolver::setMaterialPropertyVariables(){
  field->addVariable_withoutBounary("meu");
  return ;
}*/

void IdealMHDSolver::setInitialDensity(function<double(double, double)> Rho) {
    field->initializeVariable(D, Rho);
    return ;
}

void IdealMHDSolver::setInitialPressure(function<double(double, double)> Pr) {
    field->initializeVariable(P, Pr);
    return ;
}

void IdealMHDSolver::setInitialTemperature(function<double(double, double)> Tp) {
    field->initializeVariable(T, Tp);
    return ;
}

void IdealMHDSolver::setInitialVelocity(function<double(double, double)> U, function<double(double, double)> V) {
    field->initializeVariable(Vx, U);
    field->initializeVariable(Vy, V);
    return ;
}

void IdealMHDSolver::setInitialMagneticField(function<double(double, double)> BX, function<double(double, double)> BY, function<double(double, double)> BZ) {
    field->initializeVariable(Bx, BX);
    field->initializeVariable(By, BY);
    field->initializeVariable(Bz, BZ);

    return ;
}

void IdealMHDSolver::setBoundaryCondtions(string type1, string type2, string type3, string type4) {
      // Setting BC as Outflow type to test Methods
    field->setBottomBoundary(type1);
    field->setRightBoundary(type2);
    field->setTopBoundary(type3);
    field->setLeftBoundary(type4);

    return ;
}

void IdealMHDSolver::setSolver(double _CFL, double _time, int _no_of_time_steps) {
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

void IdealMHDSolver::setXMomentum() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vx, Product, DVx);
  return ;
}

void IdealMHDSolver::setYMomentum() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vy, Product, DVy);
  return ;
}

void IdealMHDSolver::setEnergy() {
  field->setFunctionsForVariables(1.0, KE, 1.0, De, Addab, DE);
  return ;
}

void IdealMHDSolver::setInternalEnergyfromPrimitive() {
  field->setFunctionsForVariables(1.0, P, 1.0, Bx, 1.0, By, 1.0, Bz, ThermoMagPressure, Pt);
  field->setFunctionsForVariables(1.0, D, 1.0, T, 1.0, Pt, IE, De);
  return ;
}

void IdealMHDSolver::setInternalEnergy() {
  field->setFunctionsForVariables(1.0, DE, 1.0, KE, Subtract, De);
  return ;
}


void IdealMHDSolver::setKineticEnergy() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vx, 1.0, Vy, KineticEnergy, KE);
  return ;
}

void IdealMHDSolver::updateVelocity() {
  field->setFunctionsForVariables(1.0, DVx, 1.0, D, Divide, Vx);
  field->setFunctionsForVariables(1.0 ,DVy, 1.0, D, Divide, Vy);
  return ;
}

void IdealMHDSolver::updateTemperature() {
  field->setFunctionsForVariables(1.0, De, 1.0, D, Temperature, T);
  return ;
}

void IdealMHDSolver::updatePressure() {
  field->setFunctionsForVariables(1.0, D, 1.0, De, 1.0, Bx, 1.0, By, 1.0, Bz, Pressure_MHDsystem, P);
  return ;
}

void IdealMHDSolver::updateThermoMagPressure() {
  field->setFunctionsForVariables(1.0, P,1.0, Bx, 1.0, By, 1.0, Bz, TotalPressure, Pt);
  return ;
}

void IdealMHDSolver::updatePrimitiveVariables() {
  updateVelocity();
  setKineticEnergy();
  setInternalEnergy();
  updateTemperature();
  updatePressure();
  updateThermoMagPressure();
  return ;
}

void IdealMHDSolver::updateConservativeVariables() {
  setXMomentum();
  setYMomentum();
  setKineticEnergy();
  setInternalEnergyfromPrimitive();
  setEnergy();
  return ;
}

void IdealMHDSolver::setInviscidFlux() {
  //cout << " Calling setInviscidFlux " << endl;

  DVxVx_plus_Pt_minus_BxBx    = field->addVariable_withBounary("quu_plus_Pt_minus_bxbx");
  DVxVy_minus_BxBy            = field->addVariable_withBounary("quv_minus_bxby");
  DVyVy_plus_Pt_minus_ByBy    = field->addVariable_withBounary("qvv_plus_Pt_minus_byby");
  DE_plus_Pt_Vx_minus_BxVdotB = field->addVariable_withBounary("qE_plus_Pt_u_minusbxVdotB");
  DE_plus_Pt_Vy_minus_ByVdotB = field->addVariable_withBounary("qE_plus_Pt_v_minusbyVdotB");
  VxBx_minus_BxVx_plus_Si     = field->addVariable_withBounary("ubx_minus_bxu_plus_si");    
  VxBy_minus_BxVy             = field->addVariable_withBounary("uBy_minus_bxv");
  VyBy_minus_ByVy_plus_Si     = field->addVariable_withBounary("vby_minus_byv_plus_si");
  BzVx                        = field->addVariable_withBounary("bzu");
  BzVy                        = field->addVariable_withBounary("bzv");

  ChBx                        = field->addVariable_withBounary("Chbx");
  ChBy                        = field->addVariable_withBounary("Chby");  

  //updateInviscidFlux();
  return ;
}

void IdealMHDSolver::updateInviscidFlux() {
  field->setFunctionsForVariables(1.0, DVx, 1.0, Vx, 1.0, Pt, 1.0, Bx, MHDMomentumFluxPressure, DVxVx_plus_Pt_minus_BxBx);
  field->setFunctionsForVariables(1.0, DVy, 1.0, Vy, 1.0, Pt, 1.0, By, MomentumFluxPressure, DVyVy_plus_Pt_minus_ByBy);
  field->setFunctionsForVariables(1.0, DVx, 1.0, Vy, 1.0, Bx, 1.0, By, MHDMomentumFlux, DVxVy_minus_BxBy);
  field->setFunctionsForVariables(1.0, DE, 1.0, Pt, 1.0, Vx, 1.0, Bx, 1.0, VdotB, MHDEnergyFlux, DE_plus_Pt_Vx_minus_BxVdotB);
  field->setFunctionsForVariables(1.0, DE, 1.0, P, 1.0, Vy, 1.0, By, 1.0, VdotB, MHDEnergyFlux, DE_plus_Pt_Vy_minus_ByVdotB);
  
  field->setFunctionsForVariables(1.0, Vx, 1.0, Bx, 1.0, Si, MHDMagFieldFlux, VxBx_minus_BxVx_plus_Si);
  field->setFunctionsForVariables(1.0, Vy, 1.0, By, 1.0, Si, MHDMagFieldFlux, VyBy_minus_ByVy_plus_Si);
  field->setFunctionsForVariables(1.0, Vx, 1.0, By, 1.0, Bx, 1.0, Vy, MHDMagFieldFlux2, VxBy_minus_BxVy);
  field->setFunctionsForVariables(1.0, Bz, -1.0, Vx, MHDZMagFlux, BzVx);
  field->setFunctionsForVariables(1.0, Bz, -1.0, Vy, MHDZMagFlux, BzVy);

  field->setFunctionsForVariables(1.0, Ch, 1.0, Bx, DivergenceFlux, ChBx);
  field->setFunctionsForVariables(1.0, Ch, 1.0, By, DivergenceFlux, ChBy);
  
  return ;
}

/*void IdealMHDSolver::setViscousFlux() {
  field->addVariable_withBounary("Tauxx");
  field->addVariable_withBounary("Tauxy");
  field->addVariable_withBounary("Tauyy");
  field->addVariable_withBounary("Eviscousx");
  field->addVariable_withBounary("Eviscousy");

  updateInviscidFlux();
  return ;
}*/

/*void IdealMHDSolver::updateViscousFlux() {
  field->setFunctionsForVariables("meu", "dudx", "dvdy", NormalViscousStress, "Tauxx");
  field->setFunctionsForVariables("meu", "dvdy", "dudx", NormalViscousStress, "Tauyy");
  field->setFunctionsForVariables("meu", "dudy", "dvdx", TangentialViscousStress, "Tauxy");
  field->setFunctionsForVariables(Vx, "Tauxx", Vy, "Tauxy", EnergyViscous, "Eviscousx");
  field->setFunctionsForVariables(Vx, "Tauxy", Vy, "Tauyy", EnergyViscous, "Eviscousy");
  return ;
}*/

/*void IdealMHDSolver::updatePrimitiveGradient() {
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

void IdealMHDSolver::setAuxillaryVariables() {
 // cout << " Calling setAuxillaryVariables " << endl;
  Si = field->addVariable_withBounary("si");
  VdotB = field->addVariable_withBounary("V.B");

  K1D   = field->addVariable_withoutBounary();
  K1DVx = field->addVariable_withoutBounary();
  K1DVy = field->addVariable_withoutBounary();
  K1DE  = field->addVariable_withoutBounary();
  K1Bx  = field->addVariable_withoutBounary();
  K1By  = field->addVariable_withoutBounary();
  K1Bz  = field->addVariable_withoutBounary();
  
  K2D   = field->addVariable_withoutBounary();
  K2DVx = field->addVariable_withoutBounary();
  K2DVy = field->addVariable_withoutBounary();
  K2DE  = field->addVariable_withoutBounary();
  K2Bx  = field->addVariable_withoutBounary();
  K2By  = field->addVariable_withoutBounary();
  K2Bz  = field->addVariable_withoutBounary();
  
  K3D   = field->addVariable_withoutBounary();
  K3DVx = field->addVariable_withoutBounary();
  K3DVy = field->addVariable_withoutBounary();
  K3DE  = field->addVariable_withoutBounary();
  K3Bx  = field->addVariable_withoutBounary();
  K3By  = field->addVariable_withoutBounary();
  K3Bz  = field->addVariable_withoutBounary();
  
  dbydxD = field->addVariable_withoutBounary();
  dbydyD = field->addVariable_withoutBounary();
  dbydxDVx = field->addVariable_withoutBounary();
  dbydyDVx = field->addVariable_withoutBounary();
  dbydxDVy = field->addVariable_withoutBounary();
  dbydyDVy = field->addVariable_withoutBounary();
  dbydxDE = field->addVariable_withoutBounary();
  dbydyDE = field->addVariable_withoutBounary();
  dbydxBx = field->addVariable_withoutBounary();
  dbydyBx = field->addVariable_withoutBounary();
  dbydxBy = field->addVariable_withoutBounary();
  dbydyBy = field->addVariable_withoutBounary();
  dbydxBz = field->addVariable_withoutBounary();
  dbydyBz = field->addVariable_withoutBounary();

  DAnalytical  = field->addVariable_withBounary("qAnalytic");
  VxAnalytical = field->addVariable_withBounary("uAnalytic");
  
  ZERO = field->addVariable_withoutBounary();
  
  return ;
}

void IdealMHDSolver::setEigenValues() {
 // cout << "Calling setEigenValues " << endl;
  C         = field->addVariable_withBounary("c");
  Vx_plus_C = field->addVariable_withBounary("u_plus_c");
  Vy_plus_C = field->addVariable_withBounary("v_plus_c");// Recheck formulation of eigen value !!
  
  Ch        = field->addVariable_withBounary("Ch");

  updateEigenValues();
  return ;
}

void IdealMHDSolver::updateEigenValues() {
 // cout << " Calling updateEigenValues " << endl;
  field->setFunctionsForVariables(1.0, D, 1.0, P, MHdMaxEigenValue, C);
  field->setFunctionsForVariables(1.0, Vx, 1.0, C, ModulusAdd, Vx_plus_C);
  field->setFunctionsForVariables(1.0, Vy, 1.0, C, ModulusAdd, Vy_plus_C);
  // Add updates for Ch , how to compute Ch ??
  return ;
}

void IdealMHDSolver::RK_Step1() {
  int FluxX[] = {DVx, DVxVx_plus_P, DVxVy, DE_plus_P_Vx};
  int FluxY[] = {DVy, DVxVy, DVyVy_plus_P, DE_plus_P_Vy};
  int FluxVarx[] = {Vx_plus_C};
  int FluxVary[] = {Vy_plus_C};
  int Var[] = {D, DVx, DVy, DE};
  int DbyDx[] = {dbydxD, dbydxDVx, dbydxDVy, dbydxDE};
  int DbyDy[] = {dbydyD, dbydyDVx, dbydyDVy, dbydyDE};

  field->delByDelX(FluxX, DbyDx, Var, "rusanov", FluxVarx, 4);
  field->delByDelY(FluxY, DbyDy, Var, "rusanov", FluxVary, 4);

  field->setFunctionsForVariables(-1.0, dbydxD, -1.0, dbydyD, Addab, K1D);
  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, Addab, K1DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, Addab, K1DVy);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, Addab, K1DE);
  
  field->setFunctionsForVariables(0.5*dt, K1D, 1.0, D, Addab, D);
  field->setFunctionsForVariables(0.5*dt, K1DVx, 1.0, DVx, Addab, DVx);
  field->setFunctionsForVariables(0.5*dt, K1DVy, 1.0, DVy, Addab, DVy);
  field->setFunctionsForVariables(0.5*dt, K1DE, 1.0, DE, Addab, DE);

  return;
}

void IdealMHDSolver::RK_Step2() {
  int FluxX[] = {DVx, DVxVx_plus_P, DVxVy, DE_plus_P_Vx};
  int FluxY[] = {DVy, DVxVy, DVyVy_plus_P, DE_plus_P_Vy};
  int FluxVarx[] = {Vx_plus_C};
  int FluxVary[] = {Vy_plus_C};
  int Var[] = {D, DVx, DVy, DE};
  int DbyDx[] = {dbydxD, dbydxDVx, dbydxDVy, dbydxDE};
  int DbyDy[] = {dbydyD, dbydyDVx, dbydyDVy, dbydyDE};

  field->delByDelX(FluxX, DbyDx, Var, "rusanov", FluxVarx, 4);
  field->delByDelY(FluxY, DbyDy, Var, "rusanov", FluxVary, 4);
  
  field->setFunctionsForVariables(-1.0, dbydxD, -1.0, dbydyD, Addab, K2D);
  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, Addab, K2DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, Addab, K2DVy);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, Addab, K2DE);
  
  field->setFunctionsForVariables(-1.5*dt, K1D, 2.0*dt, K2D, 1.0, D, Addabc, D);
  field->setFunctionsForVariables(-1.5*dt, K1DVx, 2.0*dt, K2DVx, 1.0, DVx, Addabc, DVx);
  field->setFunctionsForVariables(-1.5*dt, K1DVy, 2.0*dt, K2DVy, 1.0, DVy, Addabc, DVy);
  field->setFunctionsForVariables(-1.5*dt, K1DE, 2.0*dt, K2DE, 1.0, DE, Addabc, DE);

  return;
}

void IdealMHDSolver::RK_Step3() {
  int FluxX[] = {DVx, DVxVx_plus_P, DVxVy, DE_plus_P_Vx};
  int FluxY[] = {DVy, DVxVy, DVyVy_plus_P, DE_plus_P_Vy};
  int FluxVarx[] = {Vx_plus_C};
  int FluxVary[] = {Vy_plus_C};
  int Var[] = {D, DVx, DVy, DE};
  int DbyDx[] = {dbydxD, dbydxDVx, dbydxDVy, dbydxDE};
  int DbyDy[] = {dbydyD, dbydyDVx, dbydyDVy, dbydyDE};

  field->delByDelX(FluxX, DbyDx, Var, "rusanov", FluxVarx, 4);
  field->delByDelY(FluxY, DbyDy, Var, "rusanov", FluxVary, 4);

  field->setFunctionsForVariables(-1.0, dbydxD, -1.0, dbydyD, Addab, K3D);
  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, Addab, K3DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, Addab, K3DVy);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, Addab, K3DE);
  
  
  field->setFunctionsForVariables((7.0/6.0)*dt, K1D, -(4.0/3.0)*dt, K2D, (1.0/6.0)*dt, K3D, 1.0, D, Addabcd, D);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1DVx, -(4.0/3.0)*dt, K2DVx, (1.0/6.0)*dt, K3DVx, 1.0, DVx, Addabcd, DVx);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1DVy, -(4.0/3.0)*dt, K2DVy, (1.0/6.0)*dt, K3DVy, 1.0, DVy, Addabcd, DVy);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1DE, -(4.0/3.0)*dt, K2DE, (1.0/6.0)*dt, K3DE, 1.0, DE, Addabcd, DE);

  return;
}

void IdealMHDSolver::setTimeStep() {
  field->setFunctionsForCellCenterVariablesfromDomainVariables(1.0, Vx_plus_C, 1.0, Vy_plus_C, Maximum, UMax);
  field->FindTimestep(Dt, Dx, UMax, CFL);
  dt = field->FindMindt(Dt);

  return ;
}


void IdealMHDSolver::solve() {
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

    RK_Step1();
    
    RunShockDetector();
    RunLimiter();
    RunPositivityLimiter();

    field->updateBoundary(t);
    updatePrimitiveVariables();
   
    
    // Second Step of RK3
    updateInviscidFlux();
    updateEigenValues();
    
    RK_Step2();
    
    RunShockDetector();
    RunLimiter();
    RunPositivityLimiter(); 
    
    field->updateBoundary(t);
    updatePrimitiveVariables();

   // Third (Final) Step of RK3
    updateInviscidFlux();
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

void IdealMHDSolver::plot(string filename) {
    field->writeVTK(filename);
    return ;
}

void IdealMHDSolver::SetShockDetector(string _ShockDetector) {
  ShockDetector = _ShockDetector;
  SetShockDetectorVariables();
  return ;
}

void IdealMHDSolver::SetShockDetectorVariables() {
  
  if (ShockDetector == "KXRCF") {
    CellMarker  = field->addVariable_CellCentered();
    CellMarkerG = field->addVariable_withBounary("cellmarker");

    field->scal(0.0, CellMarkerG);
    
  }

  return ;
}

void IdealMHDSolver::RunShockDetector() {
  if (ShockDetector == "KXRCF") {
    Run_KXRCF();
 }

  return ;
}

void IdealMHDSolver::Run_KXRCF() {
 
  field->ResetVariables_CellCentered(CellMarker, 1.5);
  field->ResetMap_OutFlow();

  //field->updateOutFlowBoundary(Vx, Vy);
  //field->updateCellMarker(D, CellMarker);

  return ;
}



void IdealMHDSolver::SetLimiter(string _Limiter) {
  Limiter = _Limiter;
  SetLimiterVariables();
  return ;
}

void IdealMHDSolver::SetLimiterVariables() {
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

    dPdx = field->addVariable_withBounary("dPdx");
    dPdy = field->addVariable_withBounary("dPdy");
    dPdxMoment = field->addVariable_withoutBounary();
    dPdyMoment = field->addVariable_withoutBounary();
    field->scal(0.0, dPdxMoment);
    field->scal(0.0, dPdyMoment);

  }

  return ;
}

void IdealMHDSolver::RunLimiter() {
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

void IdealMHDSolver::Run_LiliaMomentLimiter(int v) {
  field->computeMoments(v, Moment, CellMarker);
  //field->computeMoments(v, ModMoment);
  field->setFunctionsForVariables(1.0, Moment, Copy, ModMoment);
  field->limitMoments(Moment, ModMoment, CellMarker, (N+1)*(N+1)-1);
  field->convertMomentToVariable(ModMoment, v, CellMarker);

  return ;
}

void IdealMHDSolver::RunPositivityLimiter() {
  
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

void IdealMHDSolver::Run_PositivityMomentLimiter(int v, unsigned Index) {
  field->computeMoments(v, Moment, CellMarker);
  field->scal(0.0, ModMoment);
  field->limitMoments(Moment, ModMoment, CellMarker, Index);
  field->convertMomentToVariable(ModMoment, v, CellMarker);

  return ;
}


void IdealMHDSolver::FindL2Norm(function<double(double, double)> Density, function<double(double, double)> U ) {
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

void IdealMHDSolver::checkPositivity() {
  field->checkPositivity(D, CellMarker, "One");
  field->checkPositivity(P, CellMarker, "Two");

}