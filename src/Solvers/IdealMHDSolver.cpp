#include "../../includes/Solvers/IdealMHDSolver.h"
#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/MaterialProperties.h"
#include "../../includes/Utilities/MHDFunctions.h"
#include "../../includes/Utilities/MathOperators.h"
#include "../../includes/Utilities/ThermodynamicFunctions.h"
#include "../../includes/Utilities/VectorOperations.h"


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
    field->setSystem("MHD");//change later to MHD 
    Dimension = 8;

   return ;
}

void IdealMHDSolver::setPrimitiveVariables(){
  D  = field->addVariable_withBounary("q");
  Vx = field->addVariable_withBounary("u");
  Vy = field->addVariable_withBounary("v");
  Vz = field->addVariable_withBounary("w");
  P  = field->addVariable_withBounary("P");
  T  = field->addVariable_withBounary("T");
  Bx  = field->addVariable_withBounary("Bx");
  By  = field->addVariable_withBounary("By");
  Bz  = field->addVariable_withBounary("Bz");

  Pt  = field->addVariable_withBounary("Pt");
  
  return ;
}

void IdealMHDSolver::setGradientPrimitiveVariables(){
  dPdx = field->addVariable_withBounary("dPdx");
  dPdy = field->addVariable_withBounary("dPdy");
  dBxdx = field->addVariable_withBounary("dBxdx");
  dBydy = field->addVariable_withBounary("dBydy");

  return ;
}


void IdealMHDSolver::setConservativeVariables(){
  DVx = field->addVariable_withBounary("qu");
  DVy = field->addVariable_withBounary("qv");
  DVz = field->addVariable_withBounary("qw");
  DE  = field->addVariable_withBounary("qE");
  // Check definition of internal and total energy !!
  De  = field->addVariable_withBounary("qe"); // Internal Energy, added just for ease of manipulating Energy, Remove later if not required !!
  KE  = field->addVariable_withBounary("KE"); // Kinetic Energy , " "
  
  field->addConservativeVariables(D);
  field->addConservativeVariables(DVx);
  field->addConservativeVariables(DVy);
  field->addConservativeVariables(DVz);
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

void IdealMHDSolver::setInitialVelocity(function<double(double, double)> U, function<double(double, double)> V, function<double(double, double)> W ) {
    field->initializeVariable(Vx, U);
    field->initializeVariable(Vy, V);
    field->initializeVariable(Vz, W);
    return ;
}

void IdealMHDSolver::setInitialMagneticField(function<double(double, double)> BX, function<double(double, double)> BY, function<double(double, double)> BZ) {
    field->initializeVariable(Bx, BX);
    field->initializeVariable(By, BY);
    field->initializeVariable(Bz, BZ);

    return ;
}

void IdealMHDSolver::setInitialSi(function<double(double, double)> SI) {
    field->initializeVariable(Si, SI);
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



void IdealMHDSolver::setXMomentum() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vx, Product, DVx);
  return ;
}

void IdealMHDSolver::setYMomentum() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vy, Product, DVy);
  return ;
}

void IdealMHDSolver::setZMomentum() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vz, Product, DVz);
  return ;
}

void IdealMHDSolver::setEnergy() {
  field->setFunctionsForVariables(1.0, KE, 1.0, De, Addab, DE);
  return ;
}

void IdealMHDSolver::setInternalEnergyfromPrimitive() {
  field->setFunctionsForVariables(1.0, Bx, 1.0, By, 1.0, Bz, 1.0, Bx, 1.0, By, 1.0, Bz, AdotB3d, BdotB);
  field->setFunctionsForVariables(1.0, D, 1.0, P, 1.0, BdotB, MHDIE, De);
  return ;
}

void IdealMHDSolver::setInternalEnergy() {
  field->setFunctionsForVariables(1.0, DE, 1.0, KE, Subtract, De);
  return ;
}


void IdealMHDSolver::setKineticEnergy() {
  field->setFunctionsForVariables(1.0, D, 1.0, Vx, 1.0, Vy, 1.0, Vz, KineticEnergy3d, KE);
  return ;
}

void IdealMHDSolver::updateVelocity() {
  field->setFunctionsForVariables(1.0, DVx, 1.0, D, Divide, Vx);
  field->setFunctionsForVariables(1.0 ,DVy, 1.0, D, Divide, Vy);
  field->setFunctionsForVariables(1.0 ,DVz, 1.0, D, Divide, Vz);
  return ;
}

void IdealMHDSolver::updateTemperature() {
  field->setFunctionsForVariables(1.0, De, 1.0, D, Temperature, T);
  return ;
}

void IdealMHDSolver::updatePressure() {
  field->setFunctionsForVariables(1.0, D, 1.0, De, 1.0, Bx, 1.0, By, 1.0, Bz, Pressure_MHD, P);
  return ;
}

void IdealMHDSolver::updateThermoMagPressure() {
  field->setFunctionsForVariables(1.0, P,1.0, Bx, 1.0, By, 1.0, Bz, TotalPressureMHD, Pt);
  return ;
}

void IdealMHDSolver::updatePrimitiveVariables() {
  updateVelocity();
  setKineticEnergy();
  setInternalEnergy();
  updateTemperature();
  updatePressure();
  //updateThermoMagPressure();
  return ;
}

void IdealMHDSolver::updateConservativeVariables() {
  setXMomentum();
  setYMomentum();
  setZMomentum();
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
  DVxVz_minus_BxBz            = field->addVariable_withBounary("quw_minus_bxbz");
  DVyVz_minus_ByBz            = field->addVariable_withBounary("qvw_minus_bybz");
  DE_plus_Pt_Vx_minus_BxVdotB = field->addVariable_withBounary("qE_plus_Pt_u_minusbxVdotB");
  DE_plus_Pt_Vy_minus_ByVdotB = field->addVariable_withBounary("qE_plus_Pt_v_minusbyVdotB");
  VxBx_minus_BxVx     = field->addVariable_withBounary("ubx_minus_bxu");    
  VxBy_minus_BxVy             = field->addVariable_withBounary("uby_minus_bxv");
  VyBx_minus_ByVx             = field->addVariable_withBounary("vbx_minus_uby");
  VyBy_minus_ByVy     = field->addVariable_withBounary("vby_minus_byv");
  BzVx_minus_BxVz             = field->addVariable_withBounary("bzu_minus_bxw");
  BzVy_minus_ByVz             = field->addVariable_withBounary("bzv_minus_byw");

  //ChBx                        = field->addVariable_withBounary("Chbx");
  //ChBy                        = field->addVariable_withBounary("Chby");  

  //updateInviscidFlux();
  return ;
}

void IdealMHDSolver::updateInviscidFlux() {
  
  field->setFunctionsForVariables(1.0, P, 1.0, Bx, 1.0, By, 1.0, Bz, TotalPressureMHD, Pt);
  field->setFunctionsForVariables(1.0, Vx, 1.0, Vy, 1.0, Vz, 1.0, Bx, 1.0, By, 1.0, Bz, AdotB3d, VdotB);
  field->setFunctionsForVariables(1.0, DVx, 1.0, Vx, 1.0, Pt, 1.0, Bx, MHDMomentumFluxPressure, DVxVx_plus_Pt_minus_BxBx);
  field->setFunctionsForVariables(1.0, DVy, 1.0, Vy, 1.0, Pt, 1.0, By, MHDMomentumFluxPressure, DVyVy_plus_Pt_minus_ByBy);
  field->setFunctionsForVariables(1.0, DVx, 1.0, Vy, 1.0, Bx, 1.0, By, MHDMomentumFlux, DVxVy_minus_BxBy);
  field->setFunctionsForVariables(1.0, DVx, 1.0, Vz, 1.0, Bx, 1.0, Bz, MHDMomentumFlux, DVxVz_minus_BxBz);
  field->setFunctionsForVariables(1.0, DVy, 1.0, Vz, 1.0, By, 1.0, Bz, MHDMomentumFlux, DVyVz_minus_ByBz);
  field->setFunctionsForVariables(1.0, DE, 1.0, Pt, 1.0, Vx, 1.0, Bx, 1.0, VdotB, MHDEnergyFlux, DE_plus_Pt_Vx_minus_BxVdotB);
  field->setFunctionsForVariables(1.0, DE, 1.0, Pt, 1.0, Vy, 1.0, By, 1.0, VdotB, MHDEnergyFlux, DE_plus_Pt_Vy_minus_ByVdotB);
  // Check formulation
  //field->setFunctionsForVariables(1.0, Si, Copy, VxBx_minus_BxVx);
  //field->setFunctionsForVariables(1.0, Si, Copy, VyBy_minus_ByVy);
  field->scal(0.0, VxBx_minus_BxVx);
  field->scal(0.0, VyBy_minus_ByVy);
  field->setFunctionsForVariables(1.0, Vx, 1.0, By, 1.0, Bx, 1.0, Vy, MHDMagFieldFlux, VxBy_minus_BxVy);
  field->setFunctionsForVariables(-1.0, VxBy_minus_BxVy, Copy, VyBx_minus_ByVx);
  field->setFunctionsForVariables(1.0, Vx, 1.0, Bz, 1.0, Bx, 1.0, Vz, MHDMagFieldFlux, BzVx_minus_BxVz);
  field->setFunctionsForVariables(1.0, Vy, 1.0, Bz, 1.0, By, 1.0, Vz, MHDMagFieldFlux, BzVy_minus_ByVz);
  
  //field->setFunctionsForVariables(1.0, Ch, 1.0, Bx, Product, ChBx);
  //field->setFunctionsForVariables(1.0, Ch, 1.0, By, Product, ChBy);
  
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

void IdealMHDSolver::updatePrimitiveGradient() {
  field->delByDelX(Bx, dBxdx, "central");
  field->delByDelY(By, dBydy, "central");

  return;
}

void IdealMHDSolver::setSourceTerms() {
  DeldotB_Bx = field->addVariable_withBounary("DeldotB_Bx");
  DeldotB_By = field->addVariable_withBounary("DeldotB_By");
  DeldotB_Bz = field->addVariable_withBounary("DeldotB_Bz");
  DeldotB_VdotB = field->addVariable_withBounary("DeldotB_VdotB");
  DeldotB_Vx = field->addVariable_withBounary("DeldotB_u");
  DeldotB_Vy = field->addVariable_withBounary("DeldotB_v");
  DeldotB_Vz = field->addVariable_withBounary("DeldotB_w");
  
  return ;
}

void IdealMHDSolver::updateSourceTerms(){
  //field->setFunctionsForVariables(1.0, dBxdx, 1.0, dBydy, Addab, DeldotB);
  field->updateDivergenceB(Bx, By, DeldotB);
  field->setFunctionsForVariables(1.0, DeldotB, 1.0, Bx, Product, DeldotB_Bx);
  field->setFunctionsForVariables(1.0, DeldotB, 1.0, By, Product, DeldotB_By);
  field->setFunctionsForVariables(1.0, DeldotB, 1.0, Bz, Product, DeldotB_Bz);
  field->setFunctionsForVariables(1.0, DeldotB, 1.0, VdotB, Product, DeldotB_VdotB);
  field->setFunctionsForVariables(1.0, DeldotB, 1.0, Vx, Product, DeldotB_Vx);
  field->setFunctionsForVariables(1.0, DeldotB, 1.0, Vy, Product, DeldotB_Vy);
  field->setFunctionsForVariables(1.0, DeldotB, 1.0, Vz, Product, DeldotB_Vz);

  return ;
}

void IdealMHDSolver::setAuxillaryVariables() {
 // cout << " Calling setAuxillaryVariables " << endl;
  VdotB = field->addVariable_withBounary("V.B");
  BdotB = field->addVariable_withBounary("B.B");
  DeldotB = field->addVariable_withBounary("Del.B");
  BdotB_minus_BzBz = field->addVariable_withBounary("B.B_minus_BzBz");
  Entropy = field->addVariable_withBounary("Entropy");

  K1D   = field->addVariable_withoutBounary();
  K1DVx = field->addVariable_withoutBounary();
  K1DVy = field->addVariable_withoutBounary();
  K1DVz = field->addVariable_withoutBounary();
  K1DE  = field->addVariable_withoutBounary();
  K1Bx  = field->addVariable_withoutBounary();
  K1By  = field->addVariable_withoutBounary();
  K1Bz  = field->addVariable_withoutBounary();
  
  K2D   = field->addVariable_withoutBounary();
  K2DVx = field->addVariable_withoutBounary();
  K2DVy = field->addVariable_withoutBounary();
  K2DVz = field->addVariable_withoutBounary();
  K2DE  = field->addVariable_withoutBounary();
  K2Bx  = field->addVariable_withoutBounary();
  K2By  = field->addVariable_withoutBounary();
  K2Bz  = field->addVariable_withoutBounary();
  
  K3D   = field->addVariable_withoutBounary();
  K3DVx = field->addVariable_withoutBounary();
  K3DVy = field->addVariable_withoutBounary();
  K3DVz = field->addVariable_withoutBounary();
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
  dbydxDVz = field->addVariable_withoutBounary();
  dbydyDVz = field->addVariable_withoutBounary();
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
  Cx         = field->addVariable_withBounary("cx");
  Cy         = field->addVariable_withBounary("cy");
  Vx_plus_C = field->addVariable_withBounary("u_plus_c");
  Vy_plus_C = field->addVariable_withBounary("v_plus_c");// Recheck formulation of eigen value !!
  
  //Ch        = field->addVariable_withBounary("Ch");

  //updateEigenValues();
  return ;
}

void IdealMHDSolver::updateEigenValues() {
 // cout << " Calling updateEigenValues " << endl;
  field->setFunctionsForVariables(1.0, Bx, 1.0, By, 1.0, Bz, 1.0, Bx, 1.0, By, 1.0, Bz, AdotB3d, BdotB);
  field->setFunctionsForVariables(1.0, Bx, 1.0, By, 0.0, Bz, 1.0, Bx, 1.0, By, 0.0, Bz, AdotB3d, BdotB_minus_BzBz);
  field->setFunctionsForVariables(1.0, D, 1.0, P, 1.0, Bx, 1.0, BdotB, MHDMaxEigenValue, Cx);
  field->setFunctionsForVariables(1.0, D, 1.0, P, 1.0, By, 1.0, BdotB, MHDMaxEigenValue, Cy);
  field->setFunctionsForVariables(1.0, Vx, 1.0, Cx, ModulusAdd, Vx_plus_C);
  field->setFunctionsForVariables(1.0, Vy, 1.0, Cy, ModulusAdd, Vy_plus_C);
  
  return ;
}

void IdealMHDSolver::RK_Step1() {
  int FluxX[] = {DVx, DVxVx_plus_Pt_minus_BxBx, DVxVy_minus_BxBy, DVxVz_minus_BxBz, DE_plus_Pt_Vx_minus_BxVdotB, VxBx_minus_BxVx, VxBy_minus_BxVy, BzVx_minus_BxVz};
  int FluxY[] = {DVy, DVxVy_minus_BxBy, DVyVy_plus_Pt_minus_ByBy, DVyVz_minus_ByBz, DE_plus_Pt_Vy_minus_ByVdotB, VyBx_minus_ByVx, VyBy_minus_ByVy, BzVy_minus_ByVz};
  int FluxVarx[] = {Vx_plus_C};
  int FluxVary[] = {Vy_plus_C};
  int Var[] = {D, DVx, DVy, DVz, DE, Bx, By, Bz};
  int DbyDx[] = {dbydxD, dbydxDVx, dbydxDVy, dbydxDVz, dbydxDE, dbydxBx, dbydxBy, dbydxBz};
  int DbyDy[] = {dbydyD, dbydyDVx, dbydyDVy, dbydyDVz, dbydyDE, dbydyBx, dbydyBy, dbydyBz};

  
  field->delByDelX(FluxX, DbyDx, Var, "rusanov", FluxVarx, 8);
  field->delByDelY(FluxY, DbyDy, Var, "rusanov", FluxVary, 8);

  field->setFunctionsForVariables(-1.0, dbydxD, -1.0, dbydyD, Addab, K1D);
  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, -1.0, DeldotB_Bx, Addabc, K1DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, -1.0, DeldotB_By, Addabc, K1DVy);
  field->setFunctionsForVariables(-1.0, dbydxDVz, -1.0, dbydyDVz, -1.0, DeldotB_Bz, Addabc, K1DVz);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, -1.0, DeldotB_VdotB, Addabc, K1DE);
  field->setFunctionsForVariables(-1.0, dbydxBx, -1.0, dbydyBx, -1.0, DeldotB_Vx, Addabc, K1Bx);
  field->setFunctionsForVariables(-1.0, dbydxBy, -1.0, dbydyBy, -1.0, DeldotB_Vy, Addabc, K1By);
  field->setFunctionsForVariables(-1.0, dbydxBz, -1.0, dbydyBz, -1.0, DeldotB_Vz, Addabc, K1Bz);
  
  field->setFunctionsForVariables(0.5*dt, K1D, 1.0, D, Addab, D);
  field->setFunctionsForVariables(0.5*dt, K1DVx, 1.0, DVx, Addab, DVx);
  field->setFunctionsForVariables(0.5*dt, K1DVy, 1.0, DVy, Addab, DVy);
  field->setFunctionsForVariables(0.5*dt, K1DVz, 1.0, DVz, Addab, DVz);
  field->setFunctionsForVariables(0.5*dt, K1DE, 1.0, DE, Addab, DE);
  field->setFunctionsForVariables(0.5*dt, K1Bx, 1.0, Bx, Addab, Bx);
  field->setFunctionsForVariables(0.5*dt, K1By, 1.0, By, Addab, By);
  field->setFunctionsForVariables(0.5*dt, K1Bz, 1.0, Bz, Addab, Bz);

  return;
}

void IdealMHDSolver::RK_Step2() {
  int FluxX[] = {DVx, DVxVx_plus_Pt_minus_BxBx, DVxVy_minus_BxBy, DVxVz_minus_BxBz, DE_plus_Pt_Vx_minus_BxVdotB, VxBx_minus_BxVx, VxBy_minus_BxVy, BzVx_minus_BxVz};
  int FluxY[] = {DVy, DVxVy_minus_BxBy, DVyVy_plus_Pt_minus_ByBy, DVyVz_minus_ByBz, DE_plus_Pt_Vy_minus_ByVdotB, VyBx_minus_ByVx, VyBy_minus_ByVy, BzVy_minus_ByVz};
  int FluxVarx[] = {Vx_plus_C};
  int FluxVary[] = {Vy_plus_C};
  int Var[] = {D, DVx, DVy, DVz, DE, Bx, By, Bz};
  int DbyDx[] = {dbydxD, dbydxDVx, dbydxDVy, dbydxDVz, dbydxDE, dbydxBx, dbydxBy, dbydxBz};
  int DbyDy[] = {dbydyD, dbydyDVx, dbydyDVy, dbydyDVz, dbydyDE, dbydyBx, dbydyBy, dbydyBz};

  
  field->delByDelX(FluxX, DbyDx, Var, "rusanov", FluxVarx, 8);
  field->delByDelY(FluxY, DbyDy, Var, "rusanov", FluxVary, 8);

  field->setFunctionsForVariables(-1.0, dbydxD, -1.0, dbydyD, Addab, K2D);
  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, -1.0, DeldotB_Bx, Addabc, K2DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, -1.0, DeldotB_By, Addabc, K2DVy);
  field->setFunctionsForVariables(-1.0, dbydxDVz, -1.0, dbydyDVz, -1.0, DeldotB_Bz, Addabc, K2DVz);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, -1.0, DeldotB_VdotB, Addabc, K2DE);
  field->setFunctionsForVariables(-1.0, dbydxBx, -1.0, dbydyBx, -1.0, DeldotB_Vx, Addabc, K2Bx);
  field->setFunctionsForVariables(-1.0, dbydxBy, -1.0, dbydyBy, -1.0, DeldotB_Vy, Addabc, K2By);
  field->setFunctionsForVariables(-1.0, dbydxBz, -1.0, dbydyBz, -1.0, DeldotB_Vz, Addabc, K2Bz);
  
  field->setFunctionsForVariables(-1.5*dt, K1D, 2.0*dt, K2D, 1.0, D, Addabc, D);
  field->setFunctionsForVariables(-1.5*dt, K1DVx, 2.0*dt, K2DVx, 1.0, DVx, Addabc, DVx);
  field->setFunctionsForVariables(-1.5*dt, K1DVy, 2.0*dt, K2DVy, 1.0, DVy, Addabc, DVy);
  field->setFunctionsForVariables(-1.5*dt, K1DVz, 2.0*dt, K2DVz, 1.0, DVz, Addabc, DVz);
  field->setFunctionsForVariables(-1.5*dt, K1DE, 2.0*dt, K2DE, 1.0, DE, Addabc, DE);
  field->setFunctionsForVariables(-1.5*dt, K1Bx, 2.0*dt, K2Bx, 1.0, Bx, Addabc, Bx);
  field->setFunctionsForVariables(-1.5*dt, K1By, 2.0*dt, K2By, 1.0, By, Addabc, By);
  field->setFunctionsForVariables(-1.5*dt, K1Bz, 2.0*dt, K2Bz, 1.0, Bz, Addabc, Bz);

  return;
}

void IdealMHDSolver::RK_Step3() {
  int FluxX[] = {DVx, DVxVx_plus_Pt_minus_BxBx, DVxVy_minus_BxBy, DVxVz_minus_BxBz, DE_plus_Pt_Vx_minus_BxVdotB, VxBx_minus_BxVx, VxBy_minus_BxVy, BzVx_minus_BxVz};
  int FluxY[] = {DVy, DVxVy_minus_BxBy, DVyVy_plus_Pt_minus_ByBy, DVyVz_minus_ByBz, DE_plus_Pt_Vy_minus_ByVdotB, VyBx_minus_ByVx, VyBy_minus_ByVy, BzVy_minus_ByVz};
  int FluxVarx[] = {Vx_plus_C};
  int FluxVary[] = {Vy_plus_C};
  int Var[] = {D, DVx, DVy, DVz, DE, Bx, By, Bz};
  int DbyDx[] = {dbydxD, dbydxDVx, dbydxDVy, dbydxDVz, dbydxDE, dbydxBx, dbydxBy, dbydxBz};
  int DbyDy[] = {dbydyD, dbydyDVx, dbydyDVy, dbydyDVz, dbydyDE, dbydyBx, dbydyBy, dbydyBz};

  
  field->delByDelX(FluxX, DbyDx, Var, "rusanov", FluxVarx, 8);
  field->delByDelY(FluxY, DbyDy, Var, "rusanov", FluxVary, 8);

  field->setFunctionsForVariables(-1.0, dbydxD, -1.0, dbydyD, Addab, K3D);
  field->setFunctionsForVariables(-1.0, dbydxDVx, -1.0, dbydyDVx, -1.0, DeldotB_Bx, Addabc, K3DVx);
  field->setFunctionsForVariables(-1.0, dbydxDVy, -1.0, dbydyDVy, -1.0, DeldotB_By, Addabc, K3DVy);
  field->setFunctionsForVariables(-1.0, dbydxDVz, -1.0, dbydyDVz, -1.0, DeldotB_Bz, Addabc, K3DVz);
  field->setFunctionsForVariables(-1.0, dbydxDE, -1.0, dbydyDE, -1.0, DeldotB_VdotB, Addabc, K3DE);
  field->setFunctionsForVariables(-1.0, dbydxBx, -1.0, dbydyBx, -1.0, DeldotB_Vx, Addabc, K3Bx);
  field->setFunctionsForVariables(-1.0, dbydxBy, -1.0, dbydyBy, -1.0, DeldotB_Vy, Addabc, K3By);
  field->setFunctionsForVariables(-1.0, dbydxBz, -1.0, dbydyBz, -1.0, DeldotB_Vz, Addabc, K3Bz);
    
  
  field->setFunctionsForVariables((7.0/6.0)*dt, K1D, -(4.0/3.0)*dt, K2D, (1.0/6.0)*dt, K3D, 1.0, D, Addabcd, D);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1DVx, -(4.0/3.0)*dt, K2DVx, (1.0/6.0)*dt, K3DVx, 1.0, DVx, Addabcd, DVx);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1DVy, -(4.0/3.0)*dt, K2DVy, (1.0/6.0)*dt, K3DVy, 1.0, DVy, Addabcd, DVy);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1DVz, -(4.0/3.0)*dt, K2DVz, (1.0/6.0)*dt, K3DVz, 1.0, DVz, Addabcd, DVz);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1DE, -(4.0/3.0)*dt, K2DE, (1.0/6.0)*dt, K3DE, 1.0, DE, Addabcd, DE);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1Bx, -(4.0/3.0)*dt, K2Bx, (1.0/6.0)*dt, K3Bx, 1.0, Bx, Addabcd, Bx);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1By, -(4.0/3.0)*dt, K2By, (1.0/6.0)*dt, K3By, 1.0, By, Addabcd, By);
  field->setFunctionsForVariables((7.0/6.0)*dt, K1Bz, -(4.0/3.0)*dt, K2Bz, (1.0/6.0)*dt, K3Bz, 1.0, Bz, Addabcd, Bz);
  
  return;
}

void IdealMHDSolver::setTimeStep() {
  field->setFunctionsForCellCenterVariablesfromDomainVariables(1.0, Vx_plus_C, 1.0, Vy_plus_C, Maximum, UMax);
  field->FindTimestep(Dt, Dx, UMax, CFL);
  dt = field->FindMindt(Dt);
  // Add updates for Ch , how to compute Ch ??
  //double umax = field->FindMaxCellCentered(UMax);
  //field->setConstant(umax*umax, Ch);

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
    updateEigenValues();

    if ( count%no_of_time_steps == 0) {
      setTimeStep();
      cout << "Time Step : " << dt << " , Time : " << t << "\n"; 
      count = 0;
    }
    updateInviscidFlux();
    //updatePrimitiveGradient();
    updateSourceTerms();

    RK_Step1();
    
    RunShockDetector();
    RunLimiter();
    RunPositivityLimiter();

    field->updateBoundary(t);
    updatePrimitiveVariables();
   
  
    // Second Step of RK3
    updateInviscidFlux();
    updateEigenValues();
    //updatePrimitiveGradient();
    updateSourceTerms();
    
    RK_Step2();
    
    RunShockDetector();
    RunLimiter();
    RunPositivityLimiter(); 
    
    field->updateBoundary(t);
    updatePrimitiveVariables();

   // Third (Final) Step of RK3
    updateInviscidFlux();
    updateEigenValues();
    //updatePrimitiveGradient();
    updateSourceTerms();

    RK_Step3();
    
    RunShockDetector();
    RunLimiter();
    RunPositivityLimiter(); 
   
    field->updateBoundary(t);
    updatePrimitiveVariables();
    
   t += dt; 
   count += 1;       
    }
    //updatePrimitiveGradient();
    //field->setFunctionsForVariables(1.0, dBxdx, 1.0, dBydy, Addab, DeldotB);


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
  //field->setFunctionsforDomainVariablesfromCellCenterVariables(1.0, CellMarker, SetAverage, CellMarkerG);
  
  return ;
}



void IdealMHDSolver::SetLimiter(string _Limiter) {
  Limiter = _Limiter;
  SetLimiterVariables();
  return ;
}

void IdealMHDSolver::SetLimiterVariables() {
  if (Limiter == "LiliaMoment") {
    Moment    = field->addVariable_withBounary("moment");
    field->scal(0.0, Moment);
    ModMoment = field->addVariable_withBounary("modmoment");
    field->scal(0.0, ModMoment);
    uMoment = field->addVariable_withoutBounary();
    vMoment = field->addVariable_withoutBounary();
    wMoment = field->addVariable_withoutBounary();
    qMoment = field->addVariable_withoutBounary();
    HMoment = field->addVariable_withoutBounary();
    BxMoment = field->addVariable_withoutBounary();
    ByMoment = field->addVariable_withoutBounary();
    BzMoment = field->addVariable_withoutBounary();
    uModMoment = field->addVariable_withoutBounary();
    vModMoment = field->addVariable_withoutBounary();
    wModMoment = field->addVariable_withoutBounary();
    qModMoment = field->addVariable_withoutBounary();
    HModMoment = field->addVariable_withoutBounary();
    BxModMoment = field->addVariable_withoutBounary();
    ByModMoment = field->addVariable_withoutBounary();
    BzModMoment = field->addVariable_withoutBounary();
    
    field->scal(0.0, uMoment);
    field->scal(0.0, vMoment);
    field->scal(0.0, wMoment);
    field->scal(0.0, qMoment);
    field->scal(0.0, HMoment);
    field->scal(0.0, BxMoment);
    field->scal(0.0, ByMoment);
    field->scal(0.0, BzMoment);
    field->scal(0.0, uModMoment);
    field->scal(0.0, vModMoment);
    field->scal(0.0, wModMoment);
    field->scal(0.0, qModMoment);
    field->scal(0.0, HModMoment);
    field->scal(0.0, BxModMoment);
    field->scal(0.0, ByModMoment);
    field->scal(0.0, BzModMoment);
    field->setVanderMandMatrix();
  }
  else if (Limiter == "CharacteristicLimiter") {
    // Reset for MHd System once Eigen matrices are derived !!
    field->setEigenMatrices(Dimension);
    field->setVanderMandMatrix();
    Moment    = field->addVariable_withBounary("moment");
    field->scal(0.0, Moment);
    ModMoment = field->addVariable_withBounary("modmoment");
    field->scal(0.0, ModMoment);
    // Change later by adding separate settings for positivity limiter!!
    uMoment = field->addVariable_withoutBounary();
    vMoment = field->addVariable_withoutBounary();
    wMoment = field->addVariable_withoutBounary();
    qMoment = field->addVariable_withoutBounary();
    HMoment = field->addVariable_withoutBounary();
    BxMoment = field->addVariable_withoutBounary();
    ByMoment = field->addVariable_withoutBounary();
    BzMoment = field->addVariable_withoutBounary();
    uModMoment = field->addVariable_withoutBounary();
    vModMoment = field->addVariable_withoutBounary();
    wModMoment = field->addVariable_withoutBounary();
    qModMoment = field->addVariable_withoutBounary();
    HModMoment = field->addVariable_withoutBounary();
    BxModMoment = field->addVariable_withoutBounary();
    ByModMoment = field->addVariable_withoutBounary();
    BzModMoment = field->addVariable_withoutBounary();
    
    field->scal(0.0, uMoment);
    field->scal(0.0, vMoment);
    field->scal(0.0, wMoment);
    field->scal(0.0, qMoment);
    field->scal(0.0, HMoment);
    field->scal(0.0, BxMoment);
    field->scal(0.0, ByMoment);
    field->scal(0.0, BzMoment);
    field->scal(0.0, uModMoment);
    field->scal(0.0, vModMoment);
    field->scal(0.0, wModMoment);
    field->scal(0.0, qModMoment);
    field->scal(0.0, HModMoment);
    field->scal(0.0, BxModMoment);
    field->scal(0.0, ByModMoment);
    field->scal(0.0, BzModMoment);  

    Char1 = field->addVariable_withoutBounary();
    Char2 = field->addVariable_withoutBounary();
    Char3 = field->addVariable_withoutBounary();
    Char4 = field->addVariable_withoutBounary();
    Char5 = field->addVariable_withoutBounary();
    Char6 = field->addVariable_withoutBounary();
    Char7 = field->addVariable_withoutBounary();
    Char8 = field->addVariable_withoutBounary();
    field->scal(0.0, Char1);
    field->scal(0.0, Char2);
    field->scal(0.0, Char3);
    field->scal(0.0, Char4);
    field->scal(0.0, Char5);
    field->scal(0.0, Char6);
    field->scal(0.0, Char7);
    field->scal(0.0, Char8);
    
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
    //updatePrimitiveVariables();
    int Var[] ={D, DVx, DVy, DVz, DE, Bx, By, Bz};
    int Mom[] = {qMoment, uMoment, vMoment, wMoment, HMoment, BxMoment, ByMoment, BzMoment};
    int ModMom[] = {qModMoment, uModMoment, vModMoment, wModMoment, HModMoment, BxModMoment, ByModMoment, BzModMoment}; 
    field->computeMoments(Var, Mom, CellMarker, 8);
    field->setFunctionsForVariables(1.0, qMoment, Copy, qModMoment);
    field->setFunctionsForVariables(1.0, uMoment, Copy, uModMoment);
    field->setFunctionsForVariables(1.0, vMoment, Copy, vModMoment);
    field->setFunctionsForVariables(1.0, wMoment, Copy, wModMoment);
    field->setFunctionsForVariables(1.0, HMoment, Copy, HModMoment);
    field->setFunctionsForVariables(1.0, BxMoment, Copy, BxModMoment);
    field->setFunctionsForVariables(1.0, ByMoment, Copy, ByModMoment);
    field->setFunctionsForVariables(1.0, BzMoment, Copy, BzModMoment);

    field->limitMoments(Mom, ModMom, CellMarker, (N+1)*(N+1)-1, 8);
    field->convertMomentToVariable(ModMom, Var, CellMarker, 8);
    //updateConservativeVariables();


  }

  else if(Limiter == "CharacteristicLimiter") {
    //cout << "Calling Characteristic " << endl;
    int AuxVarM[] = {qMoment, uMoment, vMoment, wMoment, HMoment, BxMoment, ByMoment, BzMoment, dPdxMoment, dPdyMoment};
    int AuxVarV[] = {qMoment, uMoment, vMoment, wMoment, HMoment, BxMoment, ByMoment, BzMoment};
    int AuxVarC[] = {Char1, Char2, Char3, Char4, Char5, Char6, Char7, Char8};
    int Var[] ={D, DVx, DVy, DVz, DE, Bx, By, Bz};
    field->computeMoments(Var, AuxVarV, CellMarker, 8);

    // Finding gradient of  Density //Pressure
    //updatePrimitiveVariables();
    //field->setFunctionsForVariables(1.0, P, 1.0, D, ThermoEntropy, Entropy);
    field->setFunctionsForVariables(1.0, D, 1.0, P, 1.0, BdotB, EntropyVar, Entropy );
    field->delByDelX(Entropy, dPdx, "central");
    field->delByDelY(Entropy, dPdy, "central");
    //field->scal(-1.0, dPdy);
    field->computeMoments(dPdx, dPdxMoment, CellMarker);
    field->computeMoments(dPdy, dPdyMoment, CellMarker);

    field->findEigenMatrices(AuxVarM, CellMarker);
    field->convertVariabletoCharacteristic(AuxVarV, AuxVarC, 0, CellMarker);
    field->limitMoments(AuxVarV, AuxVarC, CellMarker, 0);
    field->convertCharacteristictoVariable(AuxVarC, AuxVarV, 0, CellMarker);
    field->convertMomentToVariable(AuxVarV, Var, CellMarker, 8);
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
    int Var[] ={D, DVx, DVy, DVz, DE, Bx, By, Bz};
    int Mom[] = {qMoment, uMoment, vMoment, wMoment, HMoment, BxMoment, ByMoment, BzMoment};
    int ModMom[] = {qModMoment, uModMoment, vModMoment, wModMoment, HModMoment, BxModMoment, ByModMoment, BzModMoment}; 
   
    field->ResetVariables_CellCentered(CellMarker, 1.5);
    //field->resetPositivity(true);
    updatePrimitiveVariables();
    field->resetPositivity(false);
    checkPositivity();
    /*Run_PositivityMomentLimiter(DVx, N+2);
    Run_PositivityMomentLimiter(DVy, N+2);
    Run_PositivityMomentLimiter(DE, N+2);
    Run_PositivityMomentLimiter(D, N+2);
    Run_PositivityMomentLimiter(By, N+2);
    Run_PositivityMomentLimiter(Bx, N+2);
    Run_PositivityMomentLimiter(Bz, N+2);*/

    /*field->computeMoments(Var, Mom, CellMarker, 8);
    field->scal(0.0, qModMoment);
    field->scal(0.0, uModMoment);
    field->scal(0.0, vModMoment);
    field->scal(0.0, wModMoment);
    field->scal(0.0, HModMoment);
    field->scal(0.0, BxModMoment);
    field->scal(0.0, ByModMoment);
    field->scal(0.0, BzModMoment);
    field->limitMoments(Mom, ModMom, CellMarker, N+2, 8);
    field->convertMomentToVariable(ModMom, Var, CellMarker, 8);*/
    //updateConservativeVariables();


    /*updatePrimitiveVariables();
    field->resetPositivity(false);
    checkPositivity();*/
    /*Run_PositivityMomentLimiter(DVx, 0);
    Run_PositivityMomentLimiter(DVy, 0);
    Run_PositivityMomentLimiter(DE, 0);
    Run_PositivityMomentLimiter(D, 0);*/
    //Run_PositivityMomentLimiter(T, N+2); // If needed, else compute it later using q and P ..
    field->computeMoments(Var, Mom, CellMarker, 8);
    field->scal(0.0, qModMoment);
    field->scal(0.0, uModMoment);
    field->scal(0.0, vModMoment);
    field->scal(0.0, wModMoment);
    field->scal(0.0, HModMoment);
    field->scal(0.0, BxModMoment);
    field->scal(0.0, ByModMoment);
    field->scal(0.0, BzModMoment);
    field->limitMoments(Mom, ModMom, CellMarker, 0, 8);
    field->convertMomentToVariable(ModMom, Var, CellMarker, 8);
    //updateConservativeVariables();
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