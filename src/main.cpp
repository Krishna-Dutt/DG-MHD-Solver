#include "../includes/Solvers/EulerSolver.h"
#include "../includes/Solvers/AdvectionSolver.h"
#include <iostream>
#include <cmath>
#include <ctime>

#define R 285.0


using namespace std;


double U(double x, double y) {
  if ( y > 1e-6) return 0.025;
  return 0.0;
}

double V(double x, double y) {
  return 0.0;
}


double IDensity(double x, double y) {
  return 1.4;
}

double IPressure(double x, double y) {
  return 1.0;
}

double StateEq(double D, double T) {
  return D*R*T;
}

double ITemperature(double x, double y) {
  if ( x <= 0.5 && y <= 0.5 ) {
    return 1.1/(R*1.1); //1.0/(R*1.0) ;
  }
  else if ( x <= 0.5 && y > 0.5) {
    return 0.35/(R*0.5065);//0.4/(R*0.5197) ;
  }
  else if ( x > 0.5 && y > 0.5) {
    return 1.1/(R*1.1);//1.0/(R*1.0) ;
  }
  else {
    return 0.35/(R*0.5065);//0.4/(R*0.5197) ;
  }
}



double IE(double D, double T, double P) {
  return P/(gamma - 1.0) ;
}

double T(double IE, double D) {
  return IE*(gamma -1.0)/(D*R) ;
}

double Sound( double D, double P) {
  return sqrt(gamma*P/D) ;
}

double Pressures(double D, double IE) {
    return IE*(gamma - 1.0) ;
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
    int time_steps = 100;
    double CFL = 0.24;
    double time = 8.0;
    EulerSolver* a;
    a = new EulerSolver(50, 40, 1);
    a->setDomain(0.0, 0.0, 0.2, 0.1);
    a->setPrimitiveVariables();
    a->setConservativeVariables();
    a->setGradientPrimitiveVariables();
    a->setMaterialPropertyVariables();

    a->setInitialVelocity(U, V);
    a->setInitialDensity(IDensity);
    a->setInitialPressure(IPressure);
    a->setInitialTemperature(ITemperature);
    a->updateConservativeVariables(IE);

    a->setBoundaryCondtions("noslipWall", "outflow", "neumann", "inflow");
    //a->SetShockDetector("KXRCF");
    //a->SetLimiter("LiliaMoment");
    a->setSolver(CFL, time, time_steps);
    a->solve( Sound,T, Pressures, IE);
    a->FindL2Norm(IDensity, U);
    a->plot("ViscousBL_N1_50x40_Re9050_t8.vtk");
    

    delete a;
    cout << "Time Taken :: "<< (double)(clock() - tstart)/CLOCKS_PER_SEC <<"\n";
    return 0 ;
}
