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

void EulerSolver::setConservativeVariables(){
  field->addVariable_withBounary("qu");
  field->addVariable_withBounary("qv");
  field->addVariable_withBounary("qE");
  field->addVariable_withBounary("qe"); // Internal Energy, added just for ease of manipulating Energy, Remove later if not required !!
  field->addVariable_withBounary("KE"); // Kinetic Energy , " "
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

void EulerSolver::setBoundaryCondtions(string type) {
    field->setBoundaryConditions(type);
    // Setting BC as Outflow type to test Methods
    /*setBoundary("q", "periodic", "neumann", "periodic", "neumann");
    setBoundary("qu", "periodic", "neumann", "periodic", "neumann");
    setBoundary("qv", "periodic", "neumann", "periodic", "neumann");
    setBoundary("qE", "periodic", "neumann", "periodic", "neumann");*/
    return ;
}

void EulerSolver::setSolver(double _dt, double _no_of_time_steps) {
    dt = _dt;
    no_of_time_steps = _no_of_time_steps;
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
  // Mass flux qu and qv already accounted for in the conservative Variables.
  // Assuming, Requiring all Conservative and Primitive Variables have been updated !!
 // cout << " Calling updateInviscidFlux " << endl;
  field->setFunctionsForVariables("qu", "u", "P", MomentumFluxPressure, "quu_plus_P");
  field->setFunctionsForVariables("qv", "v", "P", MomentumFluxPressure, "qvv_plus_P");
  field->setFunctionsForVariables("qu", "v", MomentumFlux, "quv");
  field->setFunctionsForVariables("qE", "P", "u", EnergyFlux, "qE_plus_P_u");
  field->setFunctionsForVariables("qE", "P", "v", EnergyFlux, "qE_plus_P_v");
  return ;
}

void EulerSolver::setAuxillaryVariables() {
 // cout << " Calling setAuxillaryVariables " << endl;
  field->addVariable_withoutBounary("k1q");
  field->addVariable_withoutBounary("dqudx");
  field->addVariable_withoutBounary("dqvdy");
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
  return ;
}

void EulerSolver::setEigenValues(function<double(double,double)> SoundSpeed) {
 // cout << "Calling setEigenValues " << endl;
  field->addVariable_onlyBounary("c");
  field->addVariable_onlyBounary("u_plus_c");
  field->addVariable_onlyBounary("v_plus_c");// Recheck formulation of eigen value !!

  updateEigenValues(SoundSpeed);
  return ;
}

void EulerSolver::updateEigenValues(function<double(double,double)> SoundSpeed) {
 // cout << " Calling updateEigenValues " << endl;
  field->setFunctionsForBoundaryVariables("q", "P", SoundSpeed, "c");
  field->setFunctionsForBoundaryVariables("u", "c", ModulusAdd, "u_plus_c");
  field->setFunctionsForBoundaryVariables("v", "c", ModulusAdd, "v_plus_c");
  return ;
}

void EulerSolver::RK_Step1(string Var, string FluxX, string FluxY, string K) {
  field->delByDelX(FluxX, "dqudx", Var, "rusanov", "u_plus_c");
  field->delByDelY(FluxY, "dqvdy", Var, "rusanov", "v_plus_c");

  field->scal(0.0, K);
  field->axpy(-1.0, "dqudx", K);
  field->axpy(-1.0, "dqvdy", K);
  field->updateBoundaryVariables(Var);
  
  field->axpy(0.5*dt, K, Var);
  return;
}

void EulerSolver::RK_Step2(string Var, string FluxX, string FluxY, string K1, string K2) {
  field->delByDelX(FluxX, "dqudx", Var, "rusanov", "u_plus_c");
  field->delByDelY(FluxY, "dqvdy", Var, "rusanov", "v_plus_c");

  field->scal(0.0, K2);
  field->axpy(-1.0, "dqudx", K2);
  field->axpy(-1.0, "dqvdy", K2);
  field->updateBoundaryVariables(Var);
  
  field->axpy(-1.5*dt, K1, Var);
  field->axpy(2.0*dt, K2, Var);
  return;
}

void EulerSolver::RK_Step3(string Var, string FluxX, string FluxY, string K1, string K2, string K3) {
  field->delByDelX(FluxX, "dqudx", Var, "rusanov", "u_plus_c");
  field->delByDelY(FluxY, "dqvdy", Var, "rusanov", "v_plus_c");

  field->scal(0.0, K3);
  field->axpy(-1.0, "dqudx", K3);
  field->axpy(-1.0, "dqvdy", K3);
  field->updateBoundaryVariables(Var);
  
  field->axpy((7.0/6.0)*dt, K1, Var);
  field->axpy(-(4.0/3.0)*dt, K2, Var);
  field->axpy((1.0/6.0)*dt, K3, Var);
  return;
}


void EulerSolver::solve(function<double(double,double)> SoundSpeed ,function<double(double,double)> T, function<double(double,double)> P, function<double(double,double,double)> IE) {
  // Requires all Primitive and Conservative Variables to be setup and initialised.
  setAuxillaryVariables();
  setInviscidFlux();
  setEigenValues(SoundSpeed);

  // Till Now all variables have to be initialised !!
  // For loop to march in time !!
  for(int i=0; i < no_of_time_steps; i++) {
    // First Step of RK3
    updateInviscidFlux();
    updateEigenValues(SoundSpeed);
    
    // Mass
    RK_Step1("q", "qu", "qv", "k1q");
    // X Momentum
    RK_Step1("qu","quu_plus_P","quv","k1qu");
    // Y Momentum
    RK_Step1("qv","quv","qvv_plus_P","k1qv");
    // Energy
    RK_Step1("qE","qE_plus_P_u","qE_plus_P_v","k1qE");
    
    RunShockDetector();
    RunLimiter();
    updatePrimitiveVariables(T, P);
    /*RunShockDetector();
    RunLimiter();
    updateConservativeVariables(IE);*/
    
    
    // Second Step of RK3
    updateInviscidFlux();
    updateEigenValues(SoundSpeed);

    // Mass
    RK_Step2("q", "qu", "qv", "k1q", "k2q");
    // X Momentum
    RK_Step2("qu","quu_plus_P","quv","k1qu", "k2qu");
    // Y Momentum
    RK_Step2("qv","quv","qvv_plus_P","k1qv", "k2qv");
    // Energy
    RK_Step2("qE","qE_plus_P_u","qE_plus_P_v","k1qE", "k2qE");
    
    RunShockDetector();
    RunLimiter();
    updatePrimitiveVariables(T, P);
    RunShockDetector();
    RunLimiter();
    updateConservativeVariables(IE);
    
    

   // Third (Final) Step of RK3
    updateInviscidFlux();
    updateEigenValues(SoundSpeed);

    // Mass
    RK_Step3("q", "qu", "qv", "k1q", "k2q", "k3q");
    // X Momentum
    RK_Step3("qu","quu_plus_P","quv","k1qu", "k2qu", "k3qu");
    // Y Momentum
    RK_Step3("qv","quv","qvv_plus_P","k1qv", "k2qv", "k3qv");
    // Energy
    RK_Step3("qE","qE_plus_P_u","qE_plus_P_v","k1qE", "k2qE", "k3qE");
 
    
   RunShockDetector();
   RunLimiter();
   updatePrimitiveVariables(T, P); 
   /*RunShockDetector();
   RunLimiter();
   updateConservativeVariables(IE);*/
    
    
    updateInviscidFlux();

   /// RK3 is done, incrementing the time step. 
   // Updating the Primitive Variables !!
   // updatePrimitiveVariables(T, P);
   time += dt;        
    }
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
    /*
    field->addVariable_CellCentered("VariableMax");
    field->addVariable_CellCentered("OutFlowSize");
    field->addVariable_CellCentered("OutFlowFlux");
    field->addVariable_CellCentered("Radius");
    */
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
  /*
  field->ResetVariables_CellCentered("VariableMax");
  field->ResetVariables_CellCentered("OutFlowSize");
  field->ResetVariables_CellCentered("OutFlowFlux");
  field->ResetVariables_CellCentered("Radius");
  */
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
    field->addVariable_withBounary("Moment");
    field->addVariable_withBounary("ModifiedMoment");
    field->setVanderMandMatrix();
  }

  return ;
}

void EulerSolver::RunLimiter() {
  if ( Limiter == "LiliaMoment") {
    Run_LiliaMomentLimiter("q");
    Run_LiliaMomentLimiter("qv");
    Run_LiliaMomentLimiter("qu");
    Run_LiliaMomentLimiter("qE");
    //Run_LiliaMomentLimiter("T"); // If needed, else compute it later using q and P ..
  }

  return ;
}

void EulerSolver::Run_LiliaMomentLimiter(string v) {
  field->computeMoments(v, "Moment");
  field->computeMoments(v, "ModifiedMoment");
  field->limitMoments("Moment", "ModifiedMoment", "CellMarker");
  field->convertMomentToVariable("ModifiedMoment", v, "CellMarker");

  return ;
}

void EulerSolver::setBoundary(string v, string Bottom, string Right, string Top, string Left) {
  field->setTopBoundary(v, Top);
  field->setBottomBoundary(v, Bottom);
  field->setRightBoundary(v, Right);
  field->setLeftBoundary(v, Left);

  field->setBoundaryNeighbours(v);
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