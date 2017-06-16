#include "../includes/Solvers/EulerSolver.h"
#include <iostream>
#include <cmath>

#define R 285.0
#define gamma 1.4

using namespace std;


double U(double x, double y) {
    return 0.0;
}

double V(double x, double y) {
    return 0.0;
}

double initial(double x, double y) {
    if (x*x + y*y <= 0.25) {
      return 1.0 ;
    }
    else {
      return 0 ;
    }
    //return (exp(-(x*x +  y*y)*16.0));
}

double IDensity(double x, double y) {
  if ( x <= 0) {
    return  1.0 ;
  }
  else {
    return 0.125 ;
  }
}

double IPressure(double x, double y) {
  if ( x <= 0) {
    return  1.0 ;
  }
  else {
    return 0.1 ;
  }
}

double StateEq(double D, double T) {
  return D*R*T;
}

double ITemperature(double x, double y) {
  if ( x <= 0) {
    return  1.0/(R*1.0) ;
  }
  else {
    return 0.1/(R*0.125) ;
  }
}

double IE(double D, double T, double P) {
  return D*T*R/(gamma - 1.0) ;
}

double T(double IE, double D) {
  return IE*(gamma -1.0)/(D*R) ;
}

double Sound( double D, double T) {
  return sqrt(D*R*T);
}


int main() {
    double dt = 1e-3;
    int time_steps = 200;
    EulerSolver* a;
    a = new EulerSolver(10, 10, 2);
    a->setDomain(-1.0, -1.0, 1.0, 1.0);

    a->setInitialVelocity(U, V);
    //a->setInitialConditions(initial);
    a->setInitialDensity(IDensity);
    a->setInitialPressure(IPressure);
    a->setInitialTemperature(ITemperature);
    a->updateConservativeVariables(IE);

    a->setBoundaryCondtions("Neumann");
    a->setSolver(dt, time_steps);
    a->solve( Sound, T, StateEq, IE);
    a->plot("output.vtk");
    
    delete a;
}
