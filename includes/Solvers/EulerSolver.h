#ifndef EULERSOLVER_H
#define EULERSOLVER_H

# include "../DG_Field_2d/DG_Field_2d.h" 
#include "../Utilities/HeaderFiles.h"

using namespace std;

class EulerSolver {
private:
    int ne_x, ne_y, N;
    double x1, y1, x2, y2;
    DG_Field_2d* field;
    double time;
    double CFL;
    double dt;
    double dx;
    int no_of_time_steps;

    string ShockDetector;
    string Limiter;

public:
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This is the class constructor. The main function is to initialize the clas number of elements and the
     * order of interpolation.
     *
     * @Param _ne_x The number of elements in the x-direction.
     * @Param _ne_y The number of elements in the y-direction.
     * @Param _N    The order of interpolation used for getting the results.
     */
    /* ----------------------------------------------------------------------------*/
    EulerSolver(int _ne_x, int _ne_y, int _N);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This is the class destructor.
     *
    */
    /* ----------------------------------------------------------------------------*/
    ~EulerSolver();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This is the function for setting the domain of the problem.
     *
     * @Param _x1 The x-coordinate of the lower left corner of the domain.
     * @Param _y1 The y-coorindate of the lower left corner of the domain.
     * @Param _x2 The x-coordinate of the upper right corner of the domain.
     * @Param _y2 The y-coordinate of the upper right corner of the domain.
     */
    /* ----------------------------------------------------------------------------*/
    void setDomain(double _x1, double _y1, double _x2, double _y2);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This is the function to set the type of the boundary condition.
     *
     * @Param type This will tell the type of boundary conditions:
     *              - "periodic"  = Periodic Boundary Condition
     *              - "dirichlet" = Dirichlet Boundary Condition
     *              - "neumann"   = Neumann Boundary Condition.
     */
    /* ----------------------------------------------------------------------------*/
    void setBoundaryCondtions(string type);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to setup Primitive Variables in the domain.
     
    /* ----------------------------------------------------------------------------*/
    void setPrimitiveVariables();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to setup Conservative Variables in the domain.
     */
    /* ----------------------------------------------------------------------------*/
    void setConservativeVariables();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to initialize the density field of the domain  
     *
     * @Param functionRho This is the function used to initialize the `q` density as an input.
    */
    /* ----------------------------------------------------------------------------*/
    void setInitialDensity(function<double(double, double)> Rho);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to initialize the Pressure field of the domain  
     *
     * @Param functionP This is the function used to initialize the `P` pressure as an input.
    */
    /* ----------------------------------------------------------------------------*/
    void setInitialPressure(function<double(double, double)> P);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to initialize the Temperature field of the domain  
     *
     * @Param functionT This is the function used to initialize the `T` temperature as an input.
    */
    /* ----------------------------------------------------------------------------*/
    void setInitialTemperature(function<double(double, double)> T);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to initialize the the veclocity of the domain  
     *
     * @Param functionU This is the function used to initialize the `U` velocity as an input.
     * @Param functionV This is the function used to initialize the `V` velocity as an input.
     */
    /* ----------------------------------------------------------------------------*/
    void setInitialVelocity(function<double(double, double)> U, function<double(double, double)> V);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to set and update momentum in X-direction.
     *
    */
    /* ----------------------------------------------------------------------------*/
    void setXMomentum();
    /* ----------------------------------------------------------------------------*/
     /**
     * @Synopsis This is the function used to set and update momentum in Y-direction.
     *
    */
    /* ----------------------------------------------------------------------------*/
    void setYMomentum();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to set and update Energy ( Internal + KE ).
     *
    */
    /* ----------------------------------------------------------------------------*/
    void setEnergy();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to set and update Internal Energy.
     *
     * @Param functionIE This is the function used to compute Internal Energy, given Density, Temperature and Pressure(??).
    */
    /* ----------------------------------------------------------------------------*/
    void setInternalEnergy(function<double(double,double,double)> IE);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to set and update Internal Energy.
     *
    */
    /* ----------------------------------------------------------------------------*/
    void setInternalEnergy();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to set and update Kinetic Energy.
     *
    */
    /* ----------------------------------------------------------------------------*/
    void setKineticEnergy();
    /* ----------------------------------------------------------------------------*/
     /**
     * @Synopsis This is the function used to  update velocity Field (U,V) from Momentum.
     *
    */
    /* ----------------------------------------------------------------------------*/
    void updateVelocity();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to  update Temperature Field from Internal Energy.
     *
     * @Param functionT This is the function used to update the Temperature, given Internal Energy and Density.
    */
    /* ----------------------------------------------------------------------------*/
    void updateTemperature(function<double(double,double)> T);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to  update Pressure Field using Equation of State.
     *
     * @Param functionP This is the function used to update Pressure , given Density and Temperature, using Equation of State.
    */
    /* ----------------------------------------------------------------------------*/
    void updatePressure(function<double(double,double)> P);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to  update Primitive Variables Field using updated Conservative Variables.
     *
     * @Param functionT This is the function used to update the Temperature, given Internal Energy and Density.
     * @Param functionP This is the function used to update Pressure , given Density and Temperature, using Equation of State.
    */
    /* ----------------------------------------------------------------------------*/
    void updatePrimitiveVariables(function<double(double,double)> T, function<double(double,double)> P);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to  update Conservative Variables Field using updated Primitive Variables.
     *
     * @Param functionIE This is the function used to update Internal energy , given Density and Temperature and Pressure(??).
    */
    /* ----------------------------------------------------------------------------*/
    void updateConservativeVariables(function<double(double,double,double)> IE);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to  set up Invscid Flux Variables in the domain.
    */
    /* ----------------------------------------------------------------------------*/
    void setInviscidFlux();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to  update Invscid Flux Variables in the domain.
    */
    /* ----------------------------------------------------------------------------*/
    void updateInviscidFlux();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to  set up all auxillary Variables needed for performing RK, Time stepping.
    */
    /* ----------------------------------------------------------------------------*/
    void setAuxillaryVariables();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to  set up Eigen Values stored only at Boundaries in the domain.
     *
     * @Param functionSoundSpeed This is the function to obtain the speed of sound
    */
    /* ----------------------------------------------------------------------------*/
    void setEigenValues(function<double(double,double)> SoundSpeed);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis This is the function used to  update Eigen Values stored only at Boundaries in the domain.
     *
     * @Param functionSoundSpeed This is the function to obtain the speed of sound
    */
    /* ----------------------------------------------------------------------------*/
    void updateEigenValues(function<double(double,double)> SoundSpeed);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis   This function is used to perform First step of RK3.
     *
     * @Param Var String corresponding to Conservative Variable.
     * @Param FluxX String corresponding to flux in x direction.
     * @Param FluxY String corresponding to flux in y direction.
     * @Param K Coefficient K1
     */
    /* ----------------------------------------------------------------------------*/
    void RK_Step1(string Var, string FluxX, string FluxY, string K);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis   This function is used to perform Second step of RK3.
     *
     * @Param Var String corresponding to Conservative Variable.
     * @Param FluxX String corresponding to flux in x direction.
     * @Param FluxY String corresponding to flux in y direction.
     * @Param K1 Coefficient K1
     * @Param K2 Coefficient K2
     */
    /* ----------------------------------------------------------------------------*/
    void RK_Step2(string Var, string FluxX, string FluxY, string K1, string K2);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis   This function is used to perform Third step of RK3.
     *
     * @Param Var String corresponding to Conservative Variable.
     * @Param FluxX String corresponding to flux in x direction.
     * @Param FluxY String corresponding to flux in y direction.
     * @Param K1 Coefficient K1
     * @Param K2 Coefficient K2
     * @Param K3 Coefficient K3
     */
    /* ----------------------------------------------------------------------------*/
    void RK_Step3(string Var, string FluxX, string FluxY, string K1, string K2, string K3);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis   This function is used to set important solver parameters like dt, and no. of time steps.
     *
     * @Param _CFL CFL
     * @Param _time Evolution time.
     */
    /* ----------------------------------------------------------------------------*/
    void setSolver(double _CFL, double _time);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis   This function is find the timestep.
     *
     */
    /* ----------------------------------------------------------------------------*/
    void setTimeStep();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function does all the main functionalitites of the solver. This must be called in order to solve
     * the problem
     *
     * @Param functionSoundSpeed This is the function to obtain the speed of sound at the boundaries, given Density and Temperature.
     * @Param functionT This is the function to obtain Temperature, given Internal Energy and Density.
     * @Param functionP This is the function to obtain Pressure, given Density and Temperature, using Equation of state.
     * @Param functionIE This is the function to obtain Internal energy, given Density, Temperature and Pressure (??).
     */
    /* ----------------------------------------------------------------------------*/
    void solve(function<double(double,double)> SoundSpeed ,function<double(double,double)> T, function<double(double,double)> P, function<double(double,double,double)> IE);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function plots the function in vtk fileformat which can be further read by software packages like
     * ParaView.
     *
     * @Param filename This is the filename with which the information is to be stored.
     */
    /* ----------------------------------------------------------------------------*/
    void plot(string filename);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function sets the type of ShockDetector to be used for Euler Solver
     *
     * @Param _shockdetector This is the name of the shock detector to be used in the euler solver.
     */
    /* ----------------------------------------------------------------------------*/
    void SetShockDetector(string _shockdetector);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function sets up all additional variables required for the given Shock Detector.
     *
     */
    /* ----------------------------------------------------------------------------*/
    void SetShockDetectorVariables();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function Runs the required Shock Detector Method.
     *
     */
    /* ----------------------------------------------------------------------------*/
    void RunShockDetector();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function Resets and runs the KXRCF Shock Detector.
     *
     */
    /* ----------------------------------------------------------------------------*/
    void Run_KXRCF();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function sets the type of Limiter to be used for Euler Solver
     *
     * @Param _limiter This is the name of the limiter to be used in the euler solver.
     */
    /* ----------------------------------------------------------------------------*/
    void SetLimiter(string _limiter);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function sets up all additional variables required for the given Limiter.
     *
     */
    /* ----------------------------------------------------------------------------*/
    void SetLimiterVariables();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function Runs the required Limiter Method.
     *
     */
    /* ----------------------------------------------------------------------------*/
    void RunLimiter();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function runs the Lilia Moment Limiter.
     *
     * @Param v The variable on which the the limiter is to be performed 
     */
    /* ----------------------------------------------------------------------------*/
    void Run_LiliaMomentLimiter(string v);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function Runs the required Positivity Limiter Method.
     *
     * @Param functionT This is the function to obtain Temperature, given Internal Energy and Density.
     * @Param functionP This is the function to obtain Pressure, given Density and Temperature, using Equation of state.
     */
    /* ----------------------------------------------------------------------------*/
    void RunPositivityLimiter(function<double(double,double)> T, function<double(double,double)> P);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function runs  Positivity Limiter, using the Lilia Moment Limiter.
     *
     * @Param v The variable on which the the limiter is to be performed
     * @Param Index This is the index at which to start the limiting process. 
     */
    /* ----------------------------------------------------------------------------*/
    void Run_PositivityMomentLimiter(string v, unsigned Index);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function performs the positivity check.
     *
     */
    /* ----------------------------------------------------------------------------*/
    void checkPositivity();
    
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis   This function is set the Boundary conditions for a given variable.
     *
     * @Param Var String corresponding to the Variable.
     * @Param Bottom String corresponding to  Bottom Boundary type.
     * @Param Right String corresponding to Right Boundary type.
     * @Param Top String corresponding to Top Boundary type.
     * @Param Left String corresponding to Left Boundary type.
     */
    /* ----------------------------------------------------------------------------*/
    void setBoundary(string Var, string Bottom, string Right, string Top, string Left);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function finds the L2 norm of the converged solution ,
     *
     * @Param functionD This function gives analytical solution of Density.
     * @Param functionV This function gives analytical solution of Velocity.
     */
    /* ----------------------------------------------------------------------------*/
    void FindL2Norm(function<double(double, double)> D, function<double(double, double)> U);
 
};

#endif
