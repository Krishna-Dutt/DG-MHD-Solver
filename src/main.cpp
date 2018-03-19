#include "../includes/Solvers/IdealMHDSolver.h"
#include "../includes/Solvers/EulerSolver.h"
#include "../includes/Solvers/NavierStokesSolver.h"
#include "../includes/Solvers/AdvectionSolver.h"
#include <iostream>
#include <cmath>
#include <ctime>

#define PARALLEL false

using namespace std;

double U(double x, double y) {
  if ( x < 1e-6) return 0.135*1;//return (0.27/(0.025*0.025)) * y *( 0.05 - y);
  return 0.0;
}

double V(double x, double y) {
  return 0.0;
}

double IDensity(double x, double y) {
  return 1.0;
}

double IPressure(double x, double y) {
  if( x < 0.999) return 7.142857142857143e+01;
  return 69.70057142;
  //return 7.142857142857143e+01 + x * (69.70057142 - 7.142857142857143e+01 );
}

double StateEq(double D, double T) {
  return D*R*T;
}

double ITemperature(double x, double y) {
  return 1.0/(R*1.4) ;
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

int main(int argc, char **argv) {
    if(PARALLEL) omp_set_num_threads(8);
    clock_t tstart = clock();
    //double dt = 0.5e-3;
    int time_steps = 1;
    double CFL = 0.3;
    double time = 0.10;
    NSSolver* a;
    a = new NSSolver(80, 20, 1);
    a->setDomain(0.0, 0.0, 1.0, 0.025);
    a->setBoundaryCondtions("noslipWall", "outflow", "slipWall", "inflow");
    a->setSolver(CFL, time, time_steps);
    a->setPrimitiveVariables();
    a->setConservativeVariables();
    a->setInviscidFlux();
    a->setEigenValues();
    //a->setSourceTerms();
    a->setAuxillaryVariables();
    a->setGradientPrimitiveVariables();
    a->setViscousFlux();

    a->setInitialVelocity(U, V);
    a->setInitialDensity(IDensity);
    a->setInitialPressure(IPressure);
    a->setInitialTemperature(ITemperature);
    a->updateConservativeVariables();

    a->SetShockDetector("KXRCF");
    a->SetLimiter("LiliaMoment");
    //a->SetLimiter("CharacteristicLimiter");
    a->solve();
    a->FindL2Norm(IDensity, U);
    a->plot("PipeTest.vtk");
    

    delete a;
    cout << "Time Taken :: "<< (double)(clock() - tstart)/CLOCKS_PER_SEC <<"\n";
    return 0 ;
}
