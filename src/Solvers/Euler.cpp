#include "../../includes/Solvers/EulerSolver.h"

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

void EulerSolver::setInitialVelocity(function<double(double, double)>U, function<double(double, double)>V) {
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
    return (x*y);
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
