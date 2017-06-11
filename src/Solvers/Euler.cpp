#include "../../includes/Solvers/EulerSolver.h"
#include<iostream>

EulerSolver::EulerSolver(int _ne_x, int _ne_y, int _N) {
    ne_x = _ne_x;
    ne_y = _ne_y;
    N = _N;
    time = 0.0;
}

void EulerSolver::setDomain(double _x1, double _y1, double _x2, double _y2) {
    x1 = _x1;
    y1 = _y1;
    x2 = _x2;
    y2 = _y2;
    field = new DG_Field_2d(ne_x, ne_y, N, x1, y1, x2, y2);

    setPrimitiveVariables();
    setConservativeVariables();
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
    field->initializeVariable("p", P);
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

double MomentumFluxPressure(double f, double u, double P) {
  return (f * u + P) ;
}

double MomentumFlux(double f, double u) {
  return (f * u) ;
}

double EnergyFlux(double E, double P, double u) {
  return (E + P)*u ;
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
  field->setFunctionsForVariables("q", "T", P, "P");
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
  field->addVariable_withBounary("quu_plus_P");
  field->addVariable_withBounary("quv");
  field->addVariable_withBounary("qvv_plus_P");
  field->addVariable_withBounary("qE_plus_P_u");
  field->addVariable_withBounary("qE_plus_P_v");
  return ;
}

void EulerSolver::updateInviscidFlux() {
  // Mass flux qu and qv already accounted for in the conservative Variables.
  // Assuming, Requiring all Conservative and Primitive Variables have been updated !!
  field->setFunctionsForVariables("qu", "u", "P", MomentumFluxPressure, "quu_plus_P");
  field->setFunctionsForVariables("qv", "v", "P", MomentumFluxPressure, "qvv_plus_P");
  field->setFunctionsForVariables("qu", "v", MomentumFlux, "quv");
  field->setFunctionsForVariables("qE", "P", "u", EnergyFlux, "qE_plus_P_u");
  field->setFunctionsForVariables("qE", "P", "v", EnergyFlux, "qE_plus_P_v");
  return ;
}

void EulerSolver::setAuxillaryVariables() {
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


void EulerSolver::solve() {
    field->addVariable_withoutBounary("dqdt");
    field->addVariable_withBounary("uq");
    field->addVariable_withBounary("vq");
    field->addVariable_withoutBounary("duqdx");
    field->addVariable_withoutBounary("dvqdy");

    field->addVariable_withoutBounary("k1");
    field->addVariable_withoutBounary("k2");
    field->addVariable_withoutBounary("k3");


    // Till now the variable has been initialized.
    // This for-loop is used to march progressively in time. 
    for(int i=0; i < no_of_time_steps; i++) {
        /// First step of RK3
        field->setFunctionsForVariables("u", "q", Product, "uq");
        field->setFunctionsForVariables("v", "q", Product, "vq");
        
        field->delByDelX("uq", "duqdx", "rusanov", "u");
        field->delByDelY("vq", "dvqdy", "rusanov", "v");
        
        field->scal(0.0, "k1");
        field->axpy(-1.0, "duqdx", "k1");
        field->axpy(-1.0, "dvqdy", "k1");
        
        field->axpy(0.5*dt, "k1", "q");
        
        ///Second Step of RK3
        field->setFunctionsForVariables("u", "q", Product, "uq");
        field->setFunctionsForVariables("v", "q", Product, "vq");
        
        field->delByDelX("uq", "duqdx", "rusanov", "u");
        field->delByDelY("vq", "dvqdy", "rusanov", "v");
        
        field->scal(0.0, "k2");
        field->axpy(-1.0, "duqdx", "k2");
        field->axpy(-1.0, "dvqdy", "k2");
        
        field->axpy(-1.5*dt, "k1", "q");
        field->axpy( 2.0*dt, "k2", "q");

        /// Third(&final) step of RK3
        field->setFunctionsForVariables("u", "q", Product, "uq");
        field->setFunctionsForVariables("v", "q", Product, "vq");
        
        field->delByDelX("uq", "duqdx", "rusanov", "u");
        field->delByDelY("vq", "dvqdy", "rusanov", "v");
        
        field->scal(0.0, "k3");
        field->axpy(-1.0, "duqdx", "k3");
        field->axpy(-1.0, "dvqdy", "k3");
        
        field->axpy( (7.0/6.0)*dt, "k1", "q");
        field->axpy(-(4.0/3.0)*dt, "k2", "q");
        field->axpy( (1.0/6.0)*dt, "k3", "q");
        
        /// RK3 is done, incrementing the time step. 
        time += dt;        
    }
}

void EulerSolver::plot(string filename) {
    field->writeVTK(filename);
    return ;
}
