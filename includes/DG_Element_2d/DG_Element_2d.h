/**
 * @file DG_Element_2d.h
 * @Synopsis  This is the header file for the class which defines the DG_Element_2d class. This class contains all the
 * declarations which are necessary for defining a general 2d element. This class defines all the necessities for a
 * structured 2d DG element.
 * @author Kaushik Kulkarni
 * @version 1.0
 * @date 2017-02-18
 */



/* Copyright (c) 2016, Shivasubramanian Gopalakrishnan and Kaushik Kulkarni.
  All rights reserved.

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

     -Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
     -Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS  **AS IS** AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
     */

#ifndef DG_ELEMENT_2D_H_
#define DG_ELEMENT_2D_H_

#include "../Utilities/HeaderFiles.h"
#include "../Utilities/LobattoIntegration.h"

using namespace std;

//class DG_BoundaryElement_2d;

class DG_Element_2d {
protected:
    
    int N; /// This represents the order of polynomial interpolation which is done while solving the DG equations.
    double x_start, x_end, y_start, y_end; /// Since, this is a structured rectangular element. These variables define the position and size of the element.
    double* massMatrix;
    double* derivativeMatrix_y; /// This is the matrix which is laid out in 1-D, this would help us to find t    he $\frac{d}{dy}$ of any term.
    double* derivativeMatrix_x; /// This is the matrix which is laid out in 1-D, this would help us to find t    he $\frac{d}{dx}$ of any term. 
    double* fluxMatrix_top; /// This is the flux matrix for the top edge.
    double* fluxMatrix_right; // The Flux Matrix for the right edge.
    double* fluxMatrix_bottom; /// This would be the flux term for the the bottom edge.
    double* fluxMatrix_left; /// The Flux matrix for the left edge.
    double* inverseMassMatrix; /// The inverse of the mass matrix is also created only once in the field function. And just passed to each element.
    double* transposeDerivateMatrix_x;
    double* transposeDerivateMatrix_y;
    double* vanderMandMatrix; /// Both VanderMadn and its Inverse are created once in the DG Field and passed to each element.
    double* inverseVanderMandMatrix;
     
    DG_Element_2d* topNeighbor;
    DG_Element_2d* rightNeighbor;
    DG_Element_2d* leftNeighbor;
    DG_Element_2d* bottomNeighbor;

    string System; // System of Governing Equations being solved !!

public:

    double *X, *Y;
    double dxMin, dyMin ;
    vector<double*> variable; /// This is the map which would contain all the variable that are needed by the problem.

    vector<double**> boundaryTop;
    vector<double**> boundaryRight;
    vector<double**> boundaryLeft;
    vector<double**> boundaryBottom;
    
    vector<double**> neighboringTop;
    vector<double**> neighboringRight;
    vector<double**> neighboringBottom;
    vector<double**> neighboringLeft;

    vector<int> boundaryVariables; /// This is the variable which stores the name of all the variables whose boundary and neighboring points are stored. 
    vector<int> variableOnlyAtBoundary; // This stores all the variables which are required only at the Boundaries.

    map<string, bool> OutFlow; /// Map to flag outflow boundaries of each cell.

    bool PositivityMarker;
   

    DG_Element_2d(int _N, double x1, double y1, double x2, double y2);
    virtual ~DG_Element_2d();
    void Destroy_Matrices(); // Function to destroy all the allocated DG Matrices, once the solver is done.!!

    void setSystem(string S); 

    void addVariable_withBoundary(int v);
    void addVariable_onlyBoundary( int v);
    void addVariable_withoutBoundary(int v);
    void addVariable_CellCentered(int v);

    void initializeVariable(int v, function<double(double, double)> f);
    void setNeighboringElement(char type, DG_Element_2d* neighbor );
    void setVariableNeighbors(int v);
    void ResetMap_OutFlow();
    void updateOutFlowBoundary(int u, int v);
    double updateCellMarker(int v);

    // Function to check Positivity of data stored in cell .
    bool checkPositivity(int v, string level);
    void resetPositivity();

    // Functions to support Moment Limiters.
    void computeMoments(int v, int m);
    virtual void limitMoments(int m, int modm, unsigned Index);
    void convertMomentToVariable(int m, int v);

    // Functions to manipulate Cell Centered Variables.
    void ResetVariables_CellCentered(int v, double value = 0.0);

    // Functions for various operations on the variables.
    void delByDelX(int v, int vDash, int conserVar, string fluxType, int fluxVariable);

    void delByDelY(int v, int vDash, int conserVar, string fluxType, int fluxVariable);
    
    void delByDelX(int v, int vDash, string fluxType);

    void delByDelY(int v, int vDash, string fluxType);

    // Functions to set the various operator matrices.
    void setMassMatrix(double* m);
    void setInverseMassMatrix(double* im);
    void setderivateMatrix_x(double* d);
    void setderivateMatrix_y(double* d);
    void setTopFluxMatrix(double* f);
    void setRightFluxMatrix(double* f);
    void setLeftFluxMatrix(double* f);
    void setBottomFluxMatrix(double* f);
    void setVanderMandMatrix(double* vm);
    void setInverseVanderMandMatrix(double* ivm);
    void setTransposederivateMatrix_x(double* id);
    void setTransposederivateMatrix_y(double* id);

    // Functions to apply linear operations on the variables.
    void axpy(double a, int x, int y);
    void scal(double a, int x);
    void setFunctionsForVariables(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, unsigned, double*)>, int z); 
    void setFunctionsForVariables(double a, int w, double b, int x, double c, int y, function<void(double, double*, double, double*, double, double*, unsigned, unsigned, double*)>, int z);
    void setFunctionsForVariables(double t, int a, double u, int b, double v, int c, double x, int d, function<void(double, double*, double, double*, double, double*, double, double*, unsigned, unsigned, double*)>, int z); 
     
    // Functions to apply operations only on variables stored at Boundary
    void setFunctionsForBoundaryVariables(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, unsigned, double*)>, int z); 
    void setFunctionsForBoundaryVariables(double a, int w, double b, int x, double c, int y, function<void(double, double*, double, double*, double, double*, unsigned, unsigned, double*)>, int z); 
    
    // Fucntions to compute cell centered values
    void setFunctionsForVariablesCellCentered(int x, int y, function<double(double, double, double)>, int z); 
    

    // Functions to do various other operations on the elements.
    double l2Norm(int v1, int v2);

    // Misc. Functions on for global field data 
    double FindMax(int v);
    void FindMindx(int v);
    void FindTimestep(int dt, int dx, int U, double CFL);
    double FindMindt(int dt);

    // Virtual Function for Polymorphic behaviour
    virtual void assignBoundary( string type, char b);
    virtual void setBoundaryValue(int v, string b);
    virtual void DirichletBoundary(double *Matrix, initializer_list<int> I);
    virtual void NeumannBoundary( double *Matrix, initializer_list<int> I);
    virtual void PeriodicBoundary(double *Matrix, initializer_list<int> I);

    virtual void updateDirichlet(int v, double *Matrix);
    virtual void updateNeumann(int v, double *Matrix);
    virtual void updateBoundaryVariables(int v);

    virtual double BoundaryMinMod(int m, int Index, double Alpha, DG_Element_2d* R, DG_Element_2d* L, DG_Element_2d* T, DG_Element_2d* B);

    // Reworked Methods to update Boundary ,considering the enitre system of equation rather than individual variables

    virtual void addConservativeVariables(int v);
    virtual void addConservativeVariables(vector<int> V);

    virtual void updateBoundary(double time);
        

};
#endif
