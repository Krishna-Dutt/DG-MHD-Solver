#include "../includes/Solvers/IdealMHDSolver.h"
#include "../includes/Solvers/EulerSolver.h"
#include "../includes/Solvers/NavierStokesSolver.h"
#include "../includes/Solvers/AdvectionSolver.h"
#include <iostream>
#include <cmath>
#include <ctime>

#define PARALLEL true

using namespace std;

double U(double x, double y) {
  if ( y < 1e-10 && x >= 0.2) return 0.0;
  else if ( y <= -tan(32.6*M_PI/180.0)*(x -1.2) ) return 680.525;
  return 658.653*cos(3.1081827*M_PI/180.0);
}

double V(double x, double y) {
  if ( y <= -tan(32.6*M_PI/180.0)*(x -1.2) ) return 0.0;
  return -658.653*sin(3.1081827*M_PI/180.0);
}

double IDensity(double x, double y) {
  if ( y <= -tan(32.6*M_PI/180.0)*(x -1.2)) return 7.78e-3;
  return 8.797e-3;
}

double IPressure(double x, double y) {
  if ( y <= -tan(32.6*M_PI/180.0)*(x -1.2)) return 7.78e-3*R*288.15;
  return 8.797e-3*R*302.728;
}

double StateEq(double D, double T) {
  return D*R*T;
}

double ITemperature(double x, double y) {
  if ( y <= -tan(32.6*M_PI/180.0)*(x -1.2) ) return 288.15;
  return 302.728;
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
    if(PARALLEL) omp_set_num_threads(1);
    clock_t tstart = clock();
    //double dt = 0.5e-3;
    int time_steps = 10;
    double CFL = 0.6;
    double time = 5e-2*8e-3;
    NSSolver* a;
    a = new NSSolver(60, 35, 2);
    a->setDomain(0.0, 0.0, 1.8, 1.0);
    a->setBoundaryCondtions("AdiabaticWall", "neumann", "neumann", "dirichlet");
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
    //a->SetLimiter("LiliaMoment");
    //a->SetLimiter("CharacteristicLimiter");
    a->solve();
    a->FindL2Norm(IDensity, U);
    a->plot("ShockBLInteractionTest_N2_t4eminus2_modBC.vtk");
    
    delete a;
    cout << "Time Taken :: "<< (double)(clock() - tstart)/CLOCKS_PER_SEC <<"\n";
    return 0 ;
}
