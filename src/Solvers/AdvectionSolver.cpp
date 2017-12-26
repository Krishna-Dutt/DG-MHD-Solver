#include "../../includes/Solvers/AdvectionSolver.h"
#include <iostream>

#define D 0
#define Vx 1
#define Vy 2
#define VxD 4
#define VyD 5
#define DAnalyt 11

#define DqDt 3
#define DVxdDx 6
#define DVydDy 7
#define K1 8
#define K2 9
#define K3 10
#define ZERO 12

AdvectionSolver::AdvectionSolver(int _ne_x, int _ne_y, int _N) {
    ne_x = _ne_x;
    ne_y = _ne_y;
    N = _N;
    time = 0.0;
}

void AdvectionSolver::setDomain(double _x1, double _y1, double _x2, double _y2) {
    x1 = _x1;
    y1 = _y1;
    x2 = _x2;
    y2 = _y2;
    field = new DG_Field_2d(ne_x, ne_y, N, x1, y1, x2, y2);
    field->addVariable_withBounary(D); // q 0 
    return ;
}

void AdvectionSolver::setVelocity(function<double(double, double)>U, function<double(double, double)>V) {
    field->addVariable_withBounary(Vx); // u 1
    field->addVariable_withBounary(Vy); // v 2
    field->initializeVariable(Vx, U);
    field->initializeVariable(Vy, V);
    return ;
}

void AdvectionSolver::setInitialConditions(function<double(double, double)> I) {
    field->initializeVariable(D, I);
    return ;
}

void AdvectionSolver::setBoundaryCondtions(string type) {
    field->setBoundaryConditions(type);
    return ;
}

void AdvectionSolver::setSolver(double _dt, double _no_of_time_steps) {
    dt = _dt;
    no_of_time_steps = _no_of_time_steps;
    return ;
}


void product(double a, double* x, double b, double* y, unsigned index, unsigned N, double* z) {
    for(int i=0; i < N; ++i) {
        z[i] = x[i] * y[i];
    }
    
    return ;
}

void AxpBy(double a, double* x, double b, double* y, unsigned index, unsigned N, double* z) {
    for(int i=0; i < N; ++i) {
        z[i] = a*x[i] + b*y[i];
    }

    return;

}

void AdvectionSolver::solve() {
    field->addVariable_withoutBounary(DqDt); // dqdt 3
    field->addVariable_withBounary(VxD); // uq 4
    field->addVariable_withBounary(VyD); // vq 5

    field->addVariable_withoutBounary(DVxdDx); // duqdx 6
    field->addVariable_withoutBounary(DVydDy); // dvqdy 7

    field->addVariable_withoutBounary(K1); // k1 8 
    field->addVariable_withoutBounary(K2); // k2 9
    field->addVariable_withoutBounary(K3); // k3 10


    // Till now the variable has been initialized.
    // This for-loop is used to march progressively in time. 
    for(int i=0; i < no_of_time_steps; i++) {
       /// First step of RK3
        field->setFunctionsForVariables(1.0, Vx, 1.0, D, product, VxD);
        field->setFunctionsForVariables(1.0, Vy, 1.0, D, product, VyD);
        
        field->delByDelX(VxD, DVxdDx,D, "rusanov", Vx);
        field->delByDelY(VyD, DVydDy,D, "rusanov", Vy);
        
        field->scal(0.0, K1);
        /*field->axpy(-1.0, DVxdDx, K1);
        field->axpy(-1.0, DVydDy, K1);
        
        field->axpy(0.5*dt, K1, D);
        */
        field->setFunctionsForVariables(-1.0, DVxdDx, -1.0, DVydDy, AxpBy, K1);
        //field->setFunctionsForVariables(-1.0, DVydDy, 1.0, K1, AxpBy, K1);
        field->setFunctionsForVariables(0.5*dt, K1, 1.0, D, AxpBy, D);
        
        
        ///Second Step of RK3
        field->setFunctionsForVariables(1.0, Vx, 1.0, D, product, VxD);
        field->setFunctionsForVariables(1.0, Vy, 1.0, D, product, VyD);
        
        field->delByDelX(VxD, DVxdDx, D, "rusanov", Vx);
        field->delByDelY(VyD, DVydDy, D, "rusanov", Vy);
        
        field->scal(0.0, K2);
        /*field->axpy(-1.0, DVxdDx, K2);
        field->axpy(-1.0, DVydDy, K2);
        
        field->axpy(-1.5*dt, K1, D);
        field->axpy( 2.0*dt, K2, D);
        */
        field->setFunctionsForVariables(-1.0, DVxdDx, -1.0, DVydDy, AxpBy, K2);
        field->setFunctionsForVariables(-1.5*dt, K1, 1.0, D, AxpBy, D);
        field->setFunctionsForVariables(2.0*dt, K2, 1.0, D, AxpBy, D);
        

        /// Third(&final) step of RK3
        field->setFunctionsForVariables(1.0, Vx, 1.0, D, product, VxD);
        field->setFunctionsForVariables(1.0, Vy, 1.0, D, product, VyD);
        
        field->delByDelX(VxD, DVxdDx, D, "rusanov", Vx);
        field->delByDelY(VyD, DVydDy, D, "rusanov", Vy);
        
        field->scal(0.0, K3);
        /*field->axpy(-1.0, DVxdDx, K3);
        field->axpy(-1.0, DVydDy, K3);
        
        field->axpy( (7.0/6.0)*dt, K1, D);
        field->axpy(-(4.0/3.0)*dt, K2, D);
        field->axpy( (1.0/6.0)*dt, K3, D);
        */
        field->setFunctionsForVariables(-1.0, DVxdDx, -1.0, DVydDy, AxpBy, K3);
        field->setFunctionsForVariables((7.0/6.0)*dt, K1, 1.0, D, AxpBy, D);
        field->setFunctionsForVariables(-(4.0/3.0)*dt, K2, 1.0, D, AxpBy, D);
        field->setFunctionsForVariables((1.0/6.0)*dt, K3, 1.0, D, AxpBy, D);
        
        
        /// RK3 is done, incrementing the time step. 
        time += dt;       
    }
}

void AdvectionSolver::plot(string filename) {
    field->writeVTK(filename);
    return ;
}

void AdvectionSolver::FindL2Norm(function<double(double, double)> Density ) {
  field->addVariable_withBounary(DAnalyt); // qAnalytical 11
  
  field->addVariable_withoutBounary(ZERO); // Zero 12
  field->scal(0.0,ZERO);

  field->initializeVariable(DAnalyt, Density);

  // Computing L2Norm
  double qNorm = 0.0, qAnalytic = 0.0;
  
  qNorm = field->l2Norm(D, DAnalyt);
  qAnalytic = field->l2Norm(DAnalyt, ZERO);

  cout << "L2norms :\n  Density : " << "qNorm : " << qNorm << ", L2 Error Norm : " << qNorm/qAnalytic <<"\n";
  return ;
}
