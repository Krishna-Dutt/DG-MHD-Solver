#include "../includes/Solvers/EulerSolver.h"
#include "../includes/Solvers/AdvectionSolver.h"
#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;


/*double U(double x, double y) {
  if (x < 0.5) {
    return 1.0;
  }
  else {
    return 1.0;
  }
}

double V(double x, double y) {
  if (x < 0.5) {
    return -0.5;
  }
  else {
    return -0.5;
  }
}


double IDensity(double x, double y) {
  return 1.0 + 0.2*sin(M_PI*(x+y)) ;
}

double IPressure(double x, double y) {
  return 1.0 ;
}

double StateEq(double D, double T) {
  return D*R*T;
}

double ITemperature(double x, double y) {
   return 1.0/(R*(1.0 + 0.2*sin(M_PI*(x+y)))) ;
}

// Analytical solutions of Density  at t = 0.2 secs, for 2D Smooth Problem
double AnalyticalDensity(double x, double y) {
  double t = 4.0;
  return 1.0 + 0.2*sin(M_PI*(x+y -t*(0.5) )) ;

  
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
}*/


double U(double x, double y) {
  if ( x <= 0.5 && y <= 0.5 ) {
    return 0.8939 ;
  }
  else if ( x <= 0.5 && y > 0.5) {
    return 0.8939;
  }
  else if ( x > 0.5 && y > 0.5) {
    return 0.0 ;
  }
  else {
    return 0.0 ;
  }
}

double V(double x, double y) {
  if ( x <= 0.5 && y <= 0.5 ) {
    return 0.8939 ;
  }
  else if ( x <= 0.5 && y > 0.5) {
    return 0.0;
  }
  else if ( x > 0.5 && y > 0.5) {
    return 0.0 ;
  }
  else {
    return 0.8939 ;
  }
}


double IDensity(double x, double y) {
  if ( x <= 0.5 && y <= 0.5 ) {
    return 1.1 ;
  }
  else if ( x <= 0.5 && y > 0.5) {
    return 0.5065 ;
  }
  else if ( x > 0.5 && y > 0.5) {
    return 1.1 ;
  }
  else {
    return 0.5065 ;
  }
}

double IPressure(double x, double y) {
  if ( x <= 0.5 && y <= 0.5 ) {
    return 1.1 ;
  }
  else if ( x <= 0.5 && y > 0.5) {
    return 0.35 ;
  }
  else if ( x > 0.5 && y > 0.5) {
    return 1.1 ;
  }
  else {
    return 0.35 ;
  }
}

double StateEq(double D, double T) {
  return D*R*T;
}

double ITemperature(double x, double y) {
  if ( x <= 0.5 && y <= 0.5 ) {
    return 1.0/(R*1.0) ;
  }
  else if ( x <= 0.5 && y > 0.5) {
    return 1.0/(R*2.0) ;
  }
  else if ( x > 0.5 && y > 0.5) {
    return 1.0/(R*1.0) ;
  }
  else {
    return 1.0/(R*3.0) ;
  }
}



/*double U(double x, double y) {
  if ( x <= 0.5 && y <= 0.5 ) {
    return 0.0;//0.8939;//-0.7259 ;
  }
  else if ( x <= 0.5 && y > 0.5) {
    return 0.7276;//0.8939;//-0.7259;
  }
  else if ( x > 0.5 && y > 0.5) {
    return 0.0 ;
  }
  else {
    return 0.0 ;
  }
}

double V(double x, double y) {
  if ( x <= 0.5 && y <= 0.5 ) {
    return 0.0 ; //0.8939;//-0.7259 ;
  }
  else if ( x <= 0.5 && y > 0.5) {
    return 0.0;
  }
  else if ( x > 0.5 && y > 0.5) {
    return 0.0 ;
  }
  else {
    return 0.7276 ;//0.8939;//-0.7259 ;
  }
}


double IDensity(double x, double y) {
  if ( x <= 0.5 && y <= 0.5 ) {
    return 0.8;//1.1; //1.0 ;
  }
  else if ( x <= 0.5 && y > 0.5) {
    return 1.0;//0.5065; //0.5197;
  }
  else if ( x > 0.5 && y > 0.5) {
    return 0.5313;//1.1 ;//1.0 ;
  }
  else {
    return 1.0; //0.5065; //0.5197 ;
  }
}

double IPressure(double x, double y) {
  if ( x <= 0.5 && y <= 0.5 ) {
    return 1.0;//1.1; //1.0 ;
  }
  else if ( x <= 0.5 && y > 0.5) {
    return 1.0;//0.35; //0.4;
  }
  else if ( x > 0.5 && y > 0.5) {
    return 0.4;//1.1; //1.0 ;
  }
  else {
    return 1.0;//0.35; //0.4 ;
  }
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
*/

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
    double CFL = 0.2;
    double time = 0.25;
    EulerSolver* a;
    a = new EulerSolver(200, 200, 1);
    a->setDomain(0.0, 0.0, 1.0, 1.0);
    a->setBoundaryCondtions("neumann", "neumann", "neumann", "neumann");
    a->setSolver(CFL, time, time_steps);
    a->setPrimitiveVariables();
    a->setConservativeVariables();
    //a->setGradientPrimitiveVariables();
    //a->setMaterialPropertyVariables();

    a->setInitialVelocity(U, V);
    a->setInitialDensity(IDensity);
    a->setInitialPressure(IPressure);
    a->setInitialTemperature(ITemperature);
    a->updateConservativeVariables();

    a->setInviscidFlux();
    a->setEigenValues();
   //a->setViscousFlux();
    a->setAuxillaryVariables();

    a->SetShockDetector("KXRCF");
    //a->SetLimiter("LiliaMoment");
    a->SetLimiter("CharacteristicLimiter");
    a->solve();
    a->FindL2Norm(IDensity, U);
    a->plot("Output3.vtk");
    

    delete a;
    cout << "Time Taken :: "<< (double)(clock() - tstart)/CLOCKS_PER_SEC <<"\n";
    return 0 ;
}
