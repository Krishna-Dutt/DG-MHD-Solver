#include "../../includes/Solvers/EulerSolver.h"

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
  field->addVariable_withBounary("q");
  field->addVariable_withBounary("u");
  field->addVariable_withBounary("v");
  field->addVariable_withBounary("P");
  field->addVariable_withBounary("T");
  return ;
}

void EulerSolver::setGradientPrimitiveVariables(){
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
}


void EulerSolver::setConservativeVariables(){
  field->addVariable_withBounary("qu");
  field->addVariable_withBounary("qv");
  field->addVariable_withBounary("qE");
  field->addVariable_withBounary("qe"); // Internal Energy, added just for ease of manipulating Energy, Remove later if not required !!
  field->addVariable_withBounary("KE"); // Kinetic Energy , " "

  field->addConservativeVariables("q");
  field->addConservativeVariables("qu");
  field->addConservativeVariables("qv");
  field->addConservativeVariables("qE");
  return ;
}

void EulerSolver::setMaterialPropertyVariables(){
  field->addVariable_withoutBounary("meu");
  return ;
}

void EulerSolver::setInitialDensity(function<double(double, double)> Rho) {
    field->initializeVariable("q", Rho);
    return ;
}

void EulerSolver::setInitialPressure(function<double(double, double)> P) {
    field->initializeVariable("P", P);
    return ;
}

void EulerSolver::setInitialTemperature(function<double(double, double)> T) {
    field->initializeVariable("T", T);
    return ;
}

void EulerSolver::setInitialVelocity(function<double(double, double)> U, function<double(double, double)> V) {
    field->initializeVariable("u", U);
    field->initializeVariable("v", V);
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
   dx = field->FindMindx();
  // cout << "dx : " << dx << "\n" << "CFL : " << CFL;

    return ;
}


double Product(double x, double y) {
    return (x*y) ;
}

double KineticEnergy(double rho, double u, double v) {
  return 0.5*rho*(u*u + v*v) ;
}

double Add(double x, double y) {
  return (x + y) ;
}

double Subtract(double x, double y) {
  return (x - y) ;
}

double Subtract2(double x, double y) {
  return (x ) ;
}


double Divide(double x, double y){
  if ( y != 0.0) {
    return (x / y) ;
  }
  else {
    std::cout << "Error :: Division by Zero !! \n";
  }
  return 0 ;
}

double ModulusAdd(double x, double y) {
  return ( fabs(x) + fabs(y) ) ;
}

double MomentumFluxPressure(double f, double u, double P) {
  return (f * u + P) ;
}

double MomentumFlux(double f, double u) {
  return (f * u) ;
}

double EnergyFlux(double E, double P, double u) {
  return (E + P)*u ;
}

double PressureE(double D, double e) {
  return D*e*(1.4 -1.0) ;
}

double NormalViscousStress(double meu, double x, double y) {
  return meu * ((4.0/3.0)*x - (2.0/3.0)*y) ;
}

double TangentialViscousStress(double meu, double x, double y) {
  return meu * ( x + y) ;
}

double EnergyViscous(double u, double uTau, double v, double vTau) {
  return (u * uTau + v * vTau) ;
}

double ArtificialViscosity( double x , double y) {
  double Beta = 1.5;
  double B1 = 6.0;
  double B2 = 20.5;

 // return (Beta - 1.0)*(B1/B2)* x * (1.0/500) ;
  return 1.98e-5;
}

void EulerSolver::setXMomentum() {
  field->setFunctionsForVariables("q", "u", Product, "qu");
  return ;
}

void EulerSolver::setYMomentum() {
  field->setFunctionsForVariables("q", "v", Product, "qv");
  return ;
}

void EulerSolver::setEnergy() {
  field->setFunctionsForVariables("KE", "qe", Add, "qE");
  return ;
}

void EulerSolver::setInternalEnergy(function<double(double,double,double)> IE) {
  field->setFunctionsForVariables("q", "T", "P", IE, "qe");
  return ;
}

void EulerSolver::setInternalEnergy() {
  field->setFunctionsForVariables("qE", "KE", Subtract, "qe");
  return ;
}


void EulerSolver::setKineticEnergy() {
  field->setFunctionsForVariables("q", "u", "v", KineticEnergy, "KE");
  return ;
}

void EulerSolver::updateVelocity() {
  field->setFunctionsForVariables("qu", "q", Divide, "u");
  field->setFunctionsForVariables("qv", "q", Divide, "v");
  return ;
}

void EulerSolver::updateTemperature(function<double(double,double)> T) {
  field->setFunctionsForVariables("qe", "q", T, "T");
  return ;
}

void EulerSolver::updatePressure(function<double(double,double)> P) {
  field->setFunctionsForVariables("q", "qe", P, "P");
  return ;
}

void EulerSolver::updatePrimitiveVariables(function<double(double,double)> T, function<double(double,double)> P) {
  updateVelocity();
  setKineticEnergy();
  setInternalEnergy();
  updateTemperature(T);
  updatePressure(P);
  return ;
}

void EulerSolver::updateConservativeVariables(function<double(double,double,double)> IE) {
  setXMomentum();
  setYMomentum();
  setKineticEnergy();
  setInternalEnergy(IE);
  setEnergy();
  return ;
}

void EulerSolver::setInviscidFlux() {
  //cout << " Calling setInviscidFlux " << endl;
  field->addVariable_withBounary("quu_plus_P");
  field->addVariable_withBounary("quv");
  field->addVariable_withBounary("qvv_plus_P");
  field->addVariable_withBounary("qE_plus_P_u");
  field->addVariable_withBounary("qE_plus_P_v");

  updateInviscidFlux();
  return ;
}

void EulerSolver::updateInviscidFlux() {
  field->setFunctionsForVariables("qu", "u", "P", MomentumFluxPressure, "quu_plus_P");
  field->setFunctionsForVariables("qv", "v", "P", MomentumFluxPressure, "qvv_plus_P");
  field->setFunctionsForVariables("qu", "v", MomentumFlux, "quv");
  field->setFunctionsForVariables("qE", "P", "u", EnergyFlux, "qE_plus_P_u");
  field->setFunctionsForVariables("qE", "P", "v", EnergyFlux, "qE_plus_P_v");
  return ;
}

void EulerSolver::setViscousFlux() {
  field->addVariable_withBounary("Tauxx");
  field->addVariable_withBounary("Tauxy");
  field->addVariable_withBounary("Tauyy");
  field->addVariable_withBounary("Eviscousx");
  field->addVariable_withBounary("Eviscousy");

  updateInviscidFlux();
  return ;
}

void EulerSolver::updateViscousFlux() {
  field->setFunctionsForVariables("meu", "dudx", "dvdy", NormalViscousStress, "Tauxx");
  field->setFunctionsForVariables("meu", "dvdy", "dudx", NormalViscousStress, "Tauyy");
  field->setFunctionsForVariables("meu", "dudy", "dvdx", TangentialViscousStress, "Tauxy");
  field->setFunctionsForVariables("u", "Tauxx", "v", "Tauxy", EnergyViscous, "Eviscousx");
  field->setFunctionsForVariables("u", "Tauxy", "v", "Tauyy", EnergyViscous, "Eviscousy");
  return ;
}

void EulerSolver::updatePrimitiveGradient() {
  field->delByDelX("q", "dqdx", "central");
  field->delByDelX("u", "dudx", "central");
  field->delByDelX("v", "dvdx", "central");
  field->delByDelX("P", "dPdx", "central");
  field->delByDelX("T", "dTdx", "central");

  field->delByDelY("q", "dqdy", "central");
  field->delByDelY("u", "dudy", "central");
  field->delByDelY("v", "dvdy", "central");
  field->delByDelY("P", "dPdy", "central");
  field->delByDelY("T", "dTdy", "central");
}

void EulerSolver::setAuxillaryVariables() {
 // cout << " Calling setAuxillaryVariables " << endl;
  field->addVariable_withoutBounary("k1q");
  field->addVariable_withoutBounary("dbydx");
  field->addVariable_withoutBounary("dbydy");
  field->addVariable_withoutBounary("k1qu");
  field->addVariable_withoutBounary("k1qv");
  field->addVariable_withoutBounary("k1qE");
  field->addVariable_withoutBounary("k2q");
  field->addVariable_withoutBounary("k2qu");
  field->addVariable_withoutBounary("k2qv");
  field->addVariable_withoutBounary("k2qE");
  field->addVariable_withoutBounary("k3q");
  field->addVariable_withoutBounary("k3qu");
  field->addVariable_withoutBounary("k3qv");
  field->addVariable_withoutBounary("k3qE");

  field->addVariable_withBounary("xFlux");
  field->addVariable_withBounary("yFlux");
  return ;
}

void EulerSolver::setEigenValues(function<double(double,double)> SoundSpeed) {
 // cout << "Calling setEigenValues " << endl;
  field->addVariable_withBounary("c");
  field->addVariable_withBounary("u_plus_c");
  field->addVariable_withBounary("v_plus_c");// Recheck formulation of eigen value !!

  updateEigenValues(SoundSpeed);
  return ;
}

void EulerSolver::updateEigenValues(function<double(double,double)> SoundSpeed) {
 // cout << " Calling updateEigenValues " << endl;
  field->setFunctionsForVariables("q", "P", SoundSpeed, "c");
  field->setFunctionsForVariables("u", "c", ModulusAdd, "u_plus_c");
  field->setFunctionsForVariables("v", "c", ModulusAdd, "v_plus_c");
  return ;
}

void EulerSolver::RK_Step1(string Var, string FluxX, string FluxY, string K) {
  field->delByDelX(FluxX, "dbydx", Var, "rusanov", "u_plus_c");
  field->delByDelY(FluxY, "dbydy", Var, "rusanov", "v_plus_c");

  field->scal(0.0, K);
  field->axpy(-1.0, "dbydx", K);
  field->axpy(-1.0, "dbydy", K);
  
  field->axpy(0.5*dt, K, Var);
  return;
}

void EulerSolver::RK_Step2(string Var, string FluxX, string FluxY, string K1, string K2) {
  field->delByDelX(FluxX, "dbydx", Var, "rusanov", "u_plus_c");
  field->delByDelY(FluxY, "dbydy", Var, "rusanov", "v_plus_c");

  field->scal(0.0, K2);
  field->axpy(-1.0, "dbydx", K2);
  field->axpy(-1.0, "dbydy", K2);
  
  field->axpy(-1.5*dt, K1, Var);
  field->axpy(2.0*dt, K2, Var);
  return;
}

void EulerSolver::RK_Step3(string Var, string FluxX, string FluxY, string K1, string K2, string K3) {
  field->delByDelX(FluxX, "dbydx", Var, "rusanov", "u_plus_c");
  field->delByDelY(FluxY, "dbydy", Var, "rusanov", "v_plus_c");

  field->scal(0.0, K3);
  field->axpy(-1.0, "dbydx", K3);
  field->axpy(-1.0, "dbydy", K3);
  
  field->axpy((7.0/6.0)*dt, K1, Var);
  field->axpy(-(4.0/3.0)*dt, K2, Var);
  field->axpy((1.0/6.0)*dt, K3, Var);
  return;
}

void EulerSolver::setTimeStep() {
  //double dx = field->FindMindx();
  double speed;
  speed = max(field->FindMax("u_plus_c"),field->FindMax("v_plus_c"));
  //cout << "Max Speed  :" << field->FindMax("u_plus_c") << " : " <<field->FindMax("v_plus_c") << endl;
  //cout << "dx : " <<  dx << endl;
  dt = CFL * dx/speed;
  return ;
}


void EulerSolver::solve(function<double(double,double)> SoundSpeed ,function<double(double,double)> T, function<double(double,double)> P, function<double(double,double,double)> IE) {
  // Requires all Primitive and Conservative Variables to be setup and initialised.
  double t = 0.0;
  int count = 0;
  setAuxillaryVariables();
  setInviscidFlux();
  setEigenValues(SoundSpeed);
  setViscousFlux();

  // Till Now all variables have to be initialised !!
  // For loop to march in time !!
  while(t <= time) {
    // First Step of RK3
    
    updateInviscidFlux();
    updateEigenValues(SoundSpeed);

    if ( count%no_of_time_steps == 0) {
      setTimeStep();
      cout << "Time Step : " << dt << " , Time : " << t << "\n"; 
      count = 0;
    }

    field->setFunctionsForVariables("u_plus_c", "u_plus_c", ArtificialViscosity, "meu");
    updatePrimitiveGradient();
    updateViscousFlux();
    
    
    // Mass
    RK_Step1("q", "qu", "qv", "k1q");
    // X Momentum
    field->setFunctionsForVariables("quu_plus_P", "Tauxx", Subtract, "xFlux");
    field->setFunctionsForVariables("quv", "Tauxy", Subtract, "yFlux");
    RK_Step1("qu","xFlux","yFlux","k1qu");
    // Y Momentum
    field->setFunctionsForVariables("quv", "Tauxy", Subtract, "xFlux");
    field->setFunctionsForVariables("qvv_plus_P", "Tauyy", Subtract, "yFlux");
    RK_Step1("qv","xFlux","yFlux","k1qv");
    // Energy
    field->setFunctionsForVariables("qE_plus_P_u", "Eviscousx", Subtract, "xFlux");
    field->setFunctionsForVariables("qE_plus_P_v", "Eviscousy", Subtract, "yFlux");
    RK_Step1("qE","xFlux","yFlux","k1qE");
    
    RunShockDetector();
    RunLimiter();
    //updatePrimitiveVariables(T, P);
    RunPositivityLimiter(T, P);

    field->updateBoundary(t);
    updatePrimitiveVariables(T, P);
   
    
    // Second Step of RK3
    updateInviscidFlux();
    updateEigenValues(SoundSpeed);

    field->setFunctionsForVariables("u_plus_c", "u_plus_c", ArtificialViscosity, "meu");
    updatePrimitiveGradient();
    updateViscousFlux();

    // Mass
    RK_Step2("q", "qu", "qv", "k1q", "k2q");
    // X Momentum
    field->setFunctionsForVariables("quu_plus_P", "Tauxx", Subtract, "xFlux");
    field->setFunctionsForVariables("quv", "Tauxy", Subtract, "yFlux");
    RK_Step2("qu","xFlux","yFlux","k1qu", "k2qu");
    // Y Momentum
    field->setFunctionsForVariables("quv", "Tauxy", Subtract, "xFlux");
    field->setFunctionsForVariables("qvv_plus_P", "Tauyy", Subtract, "yFlux");
    RK_Step2("qv","xFlux","yFlux","k1qv", "k2qv");
    // Energy
    field->setFunctionsForVariables("qE_plus_P_u", "Eviscousx", Subtract, "xFlux");
    field->setFunctionsForVariables("qE_plus_P_v", "Eviscousy", Subtract, "yFlux");
    RK_Step2("qE","xFlux","yFlux","k1qE", "k2qE");
    
    RunShockDetector();
    RunLimiter();
    //updatePrimitiveVariables(T, P);
    RunPositivityLimiter(T, P); 
    
    field->updateBoundary(t);
    updatePrimitiveVariables(T, P);

   // Third (Final) Step of RK3
    updateInviscidFlux();
    updateEigenValues(SoundSpeed);

    field->setFunctionsForVariables("u_plus_c", "u_plus_c", ArtificialViscosity, "meu");
    updatePrimitiveGradient();
    updateViscousFlux();

    // Mass
    RK_Step3("q", "qu", "qv", "k1q", "k2q", "k3q");
    // X Momentum
    field->setFunctionsForVariables("quu_plus_P", "Tauxx", Subtract, "xFlux");
    field->setFunctionsForVariables("quv", "Tauxy", Subtract, "yFlux");
    RK_Step3("qu","xFlux","yFlux","k1qu", "k2qu", "k3qu");
    // Y Momentum
    field->setFunctionsForVariables("quv", "Tauxy", Subtract, "xFlux");
    field->setFunctionsForVariables("qvv_plus_P", "Tauyy", Subtract, "yFlux");
    RK_Step3("qv","xFlux","yFlux","k1qv", "k2qv", "k3qv");
    // Energy
    field->setFunctionsForVariables("qE_plus_P_u", "Eviscousx", Subtract, "xFlux");
    field->setFunctionsForVariables("qE_plus_P_v", "Eviscousy", Subtract, "yFlux");
    RK_Step3("qE","xFlux","yFlux","k1qE", "k2qE", "k3qE");
 
    
   RunShockDetector();
   RunLimiter();
   //updatePrimitiveVariables(T, P);
   RunPositivityLimiter(T, P); 
   
    field->updateBoundary(t);
    updatePrimitiveVariables(T, P);

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
    field->addVariable_CellCentered("CellMarker");
    field->addVariable_withBounary("CellMarkerGlobal");
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
 
  field->ResetVariables_CellCentered("CellMarker", 1.5);
  field->ResetMap_OutFlow();

  field->updateOutFlowBoundary("u", "v");
  field->updateCellMarker("q", "CellMarker");

  return ;
}



void EulerSolver::SetLimiter(string _Limiter) {
  Limiter = _Limiter;
  SetLimiterVariables();
  return ;
}

void EulerSolver::SetLimiterVariables() {
  if (Limiter == "LiliaMoment") {
    field->addVariable_CellCentered("CellMarker");
    field->ResetVariables_CellCentered("CellMarker", 1.5);
    field->addVariable_withBounary("CellMarkerGlobal");
    field->scal(0.0, "CellMarkerGlobal");
    field->addVariable_CellCentered("Max");
    field->addVariable_CellCentered("Min");


    field->addVariable_withBounary("Moment");
    field->addVariable_withBounary("ModifiedMoment");
    field->setVanderMandMatrix();
  }

  return ;
}

void EulerSolver::RunLimiter() {
  if ( Limiter == "LiliaMoment") {
    field->resetPositivity();
    Run_LiliaMomentLimiter("qu");
    Run_LiliaMomentLimiter("qv");
    Run_LiliaMomentLimiter("qE");
    Run_LiliaMomentLimiter("q");
    //Run_LiliaMomentLimiter("T"); // If needed, else compute it later using q and P ..
  }

  return ;
}

void EulerSolver::Run_LiliaMomentLimiter(string v) {
  field->computeMoments(v, "Moment");
  field->computeMoments(v, "ModifiedMoment");
  field->limitMoments("Moment", "ModifiedMoment", "CellMarker", (N+1)*(N+1)-1);
  field->convertMomentToVariable("ModifiedMoment", v, "CellMarker");

  return ;
}

void EulerSolver::RunPositivityLimiter(function<double(double,double)> T ,function<double(double,double)> P) {
  field->resetPositivity();
  if ( Limiter == "LiliaMoment") {
    updatePrimitiveVariables(T, P);
    checkPositivity();
    Run_PositivityMomentLimiter("qu", N+2);
    Run_PositivityMomentLimiter("qv", N+2);
    Run_PositivityMomentLimiter("qE", N+2);
    Run_PositivityMomentLimiter("q", N+2);

    updatePrimitiveVariables(T, P);
    checkPositivity();
    Run_PositivityMomentLimiter("qu", 0);
    Run_PositivityMomentLimiter("qv", 0);
    Run_PositivityMomentLimiter("qE", 0);
    Run_PositivityMomentLimiter("q", 0);
    //Run_PositivityMomentLimiter("T", N+2); // If needed, else compute it later using q and P ..
  }

  return ;
}

void EulerSolver::Run_PositivityMomentLimiter(string v, unsigned Index) {
  field->computeMoments(v, "Moment");
  field->scal(0.0, "ModifiedMoment");
  field->limitMoments("Moment", "ModifiedMoment", "CellMarker", Index);
  field->convertMomentToVariable("ModifiedMoment", v, "CellMarker");

  return ;
}


void EulerSolver::FindL2Norm(function<double(double, double)> D, function<double(double, double)> U ) {
  field->addVariable_withBounary("qAnalytical");
  field->addVariable_withBounary("uAnalytical");
  
  field->addVariable_withoutBounary("Zero");
  field->scal(0.0,"Zero");

  field->initializeVariable("qAnalytical", D);
  field->initializeVariable("uAnalytical", U);

  // Computing L2Norm

  double qNorm = 0.0, qAnalytic = 0.0;
  double uNorm = 0.0, uAnalytic = 0.0;
  
  qNorm = field->l2Norm("q", "qAnalytical");
  qAnalytic = field->l2Norm("qAnalytical", "Zero");
  uNorm = field->l2Norm("u", "uAnalytical");
  uAnalytic = field->l2Norm("uAnalytical", "Zero");

  cout << "L2norms :\n  Density : " <<qNorm <<" : L2 Norm : "<< qNorm/qAnalytic <<"\n";
  cout <<" u Velocity : " <<uNorm << " : L2 Norm :  " << uNorm/uAnalytic <<"\n";
  return ;
}

void EulerSolver::checkPositivity() {
  field->checkPositivity("q", "CellMarker", "One");
  field->checkPositivity("P", "CellMarker", "Two");

}