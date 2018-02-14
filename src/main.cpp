#include "../includes/Solvers/IdealMHDSolver.h"
#include "../includes/Solvers/EulerSolver.h"
#include "../includes/Solvers/AdvectionSolver.h"
#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;


double U(double x, double y) {
  double r0 = 0.1, r1 = 0.115, u0 = 2.0;
  double r = sqrt( pow(x -0.5,2.0) + pow(y-0.5,2.0) );
  double f = (r1 - r)/(r1 - r0);

  if ( r < r0) {
    return (u0/r0)*(-y + 0.5);
  }
  else if ( r >= r0  &&  r < r1){
    return (f*u0/r)*(-y + 0.5);
  }
  else {
    return 0.0;
  }
}

double V(double x, double y) {
  double r0 = 0.1, r1 = 0.115, u0 = 2.0;
  double r = sqrt( pow(x -0.5,2.0) + pow(y-0.5,2.0) );
  double f = (r1 - r)/(r1 - r0);

  if ( r < r0) {
    return (u0/r0)*(x - 0.5);
  }
  else if ( r >= r0  &&  r < r1){
    return (f*u0/r)*(x - 0.5);
  }
  else {
    return 0.0;
  }
}

double W(double x, double y) {
  return 0.0;
}

double BX(double x, double y) {
  return 5.0/sqrt(4.0*M_PI);;
}

double BY(double x, double y) {
  return 0.0;
}

double BZ(double x, double y) {
  return 0.0 ;
}

double SI(double x, double y) {
  return 0.0 ;
}

double IDensity(double x, double y) {
  double r0 = 0.1, r1 = 0.115, u0 = 2.0;
  double r = sqrt( pow(x -0.5,2.0) + pow(y-0.5,2.0) );
  double f = (r1 - r)/(r1 - r0);

  if ( r < r0) {
    return 10;
  }
  else if ( r >= r0  &&  r < r1){
    return 1.0 + 9.0*f;
  }
  else {
    return 1.0;
  }
}

double IPressure(double x, double y) {
  return 1.0; 
}

double StateEq(double D, double T) {
  return D*R*T;
}

double ITemperature(double x, double y) {
  return 3.0/(5.0*R);
}


// Analytical solutions of Density and Pressure at t = 0.2 secs, for 1D Sod's Shock Tube
double AnalyticalDensity(double x, double y) {
  if (0.0 <= x && x <= 0.26) {
    return 1.0 ;
  }
  else if (0.26 <= x && x <= 0.485) {
    return 1.0 + (x-0.26)*(0.42-1.0)/(0.485-0.26) ;
  }
  else if (0.485 <= x && x <= 0.682) {
    return 0.423 ;
  }
  else if (0.682 < x && x <= 0.852) {
    return 0.27 ;
  }
  else {
    return 0.125 ;
  }
}

double AnalyticalVelocity(double x, double y) {
  if (0.0 <= x && x <= 0.26) {
    return 0.0 ;
  }
  else if (0.26 <= x && x <= 0.485) {
    return (x-0.26)*(0.926)/(0.485-0.26) ;
  }
  else if (0.485 <= x & x <= 0.852) {
    return 0.926 ;
  }
  else {
    return 0.0 ;
  }
}

int main() {
    clock_t tstart = clock();
    //double dt = 0.5e-3;
    int time_steps = 10;
    double CFL = 0.1;
    double time = 0.15;
    IdealMHDSolver* a;
    a = new IdealMHDSolver(200, 200, 1);
    a->setDomain(0.0, 0.0, 1.0, 1.0);
    a->setBoundaryCondtions("periodic", "periodic", "periodic", "periodic");
    a->setSolver(CFL, time, time_steps);
    a->setPrimitiveVariables();
    a->setConservativeVariables();
    a->setInviscidFlux();
    a->setEigenValues();
    a->setSourceTerms();
    a->setAuxillaryVariables();
    a->setGradientPrimitiveVariables();
    //a->setMaterialPropertyVariables();

    a->setInitialVelocity(U, V, W);
    a->setInitialDensity(IDensity);
    a->setInitialPressure(IPressure);
    a->setInitialTemperature(ITemperature);
    a->setInitialMagneticField(BX, BY, BZ);
    //a->setInitialSi(SI);
    a->updateConservativeVariables();

    a->SetShockDetector("KXRCF");
    a->SetLimiter("LiliaMoment");
    //a->SetLimiter("CharacteristicLimiter");
    a->solve();
    a->FindL2Norm(IDensity, U);
    a->plot("RotorTest.vtk");
    

    delete a;
    cout << "Time Taken :: "<< (double)(clock() - tstart)/CLOCKS_PER_SEC <<"\n";
    return 0 ;
}
