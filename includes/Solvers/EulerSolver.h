#ifndef EULERSOLVER_H
#define EULERSOLVER_H

# include "../DG_Field_2d/DG_Field_2d.h" 
#include <functional>

using namespace std;

class EulerSolver {
private:
    int ne_x, ne_y, N;
    double x1, y1, x2, y2;
    DG_Field_2d* field;
    double time;
    double dt;
    int no_of_time_steps;

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
     * @Synopsis   This function is used to set important solver parameters like dt, and no. of time steps.
     *
     * @Param _dt The time step for each iteration.
     * @Param _no_of_time_steps The number of time steps that must be used which is also the number of time iterations that must
     * be performed
     */
    /* ----------------------------------------------------------------------------*/
    void setSolver(double _dt, double _no_of_time_steps);
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function does all the main functionalitites of the solver. This must be called in order to solve
     * the problem
     */
    /* ----------------------------------------------------------------------------*/
    void solve();
    /* ----------------------------------------------------------------------------*/
    /**
     * @Synopsis  This function plots the function in vtk fileformat which can be further read by software packages like
     * ParaView.
     *
     * @Param filename This is the filename with which the information is to be stored.
     */
    /* ----------------------------------------------------------------------------*/
    void plot(string filename);
};

#endif
