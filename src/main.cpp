#include "../includes/Solvers/EulerSolver.h"
#include "../includes/Solvers/AdvectionSolver.h"
#include <iostream>
#include <cmath>
#include <ctime>

#define R 285.0
#define gamma 1.4

using namespace std;


double U(double x, double y) {
    return 1.0;
}

double V(double x, double y) {
    return 0.0;
}

double initial(double x, double y) {
   /* if (x*x + y*y <= 0.25) {
      return 1.0 ;
    }
    else {
      return 0 ;
    }*/
    //return (exp(-(x*x)*16.0));
    return -x ;
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
  return P/(gamma - 1.0) ;
}

double T(double IE, double D) {
  return IE*(gamma -1.0)/(D*R) ;
}

double Sound( double D, double P) {
  return sqrt(gamma*P/D);
}


int main() {
    clock_t tstart = clock();
    double dt = 1e-2;
    int time_steps = 100;
    AdvectionSolver* a;
    a = new AdvectionSolver(30, 1, 2);
    a->setDomain(-1.0, -1.0, 1.0, 1.0);

    a->setVelocity(U, V);
    a->setInitialConditions(initial);
    //a->setBoundaryCondtions("periodic");
    /*a->setInitialDensity(IDensity);
    a->setInitialPressure(IPressure);
    a->setInitialTemperature(ITemperature);
    a->updateConservativeVariables(IE);

    a->setBoundaryCondtions("periodicY");
    a->SetShockDetector("KXRCF");
    a->SetLimiter("LiliaMoment");*/
    a->setSolver(dt, time_steps);
    //a->solve( Sound, T, StateEq, IE);
    a->solve();
    a->plot("output.vtk");
    
    delete a;
    cout << "Time Taken :: "<< (double)(clock() - tstart)/CLOCKS_PER_SEC <<"\n";
    return 0 ;
}
