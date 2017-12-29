/**
  * @file DG_Field_2d.h
  * @Synopsis  This is the file for the class `Field`, which the data-structure for storing DG Fields, This is the
  * header file and it only stores the declarations to all the functions and the member variables of the class.
  * @author Kaushik Kulkarni
  * @version 1.0
  * @date 2017-02-18
  */

#ifndef DG_FIELD_2D_H_
#define DG_FIELD_2D_H_

#include "../Utilities/HeaderFiles.h"

using namespace std; 

class DG_Element_2d; 

class DG_Field_2d {
private:
    int N; /// The order of the polynomial of interpolation.
    int ne_x, ne_y; /// Since this is a structured grid we can define the grid with the help of number of elements in the x-direction and in the y-direction.
    double x1, y1, x2, y2;


public:
    vector< vector<DG_Element_2d*> > elements;
    vector<double*> domainVariable; // To store variables that have values at all nodes
    vector<double*> cellcenterVariable; // To store variables with only one value per element
    vector<bool> PositivityMarker;

    map<int, string> boundaryVariableNames;
    vector<int> variableNames; // This is stores all the variables which have been added to the field.
    vector<int> variablesWithBoundaryInfo; // This stores all the variables whose boundary info. is also known.
    vector<int> variableOnlyAtBoundary; // This stores all the variables who are required only at the Boundaries.

    DG_Field_2d(int _nex, int _ney, int _N, double _x1, double _y1, double _x2, double _y2);
    ~DG_Field_2d();

    void setSystem(string S);

    int addVariable_withBounary(string V);
    int addVariable_withoutBounary();
    int addVariable_onlyBounary();
    int addVariable_CellCentered();
    void initializeVariable(int v, function<double(double, double)>);

    // Functions to set Boundary Conditions
    void setBoundaryConditions(string type);
    void setBoundaryConditions(int v, string type, string b);
    void setTopBoundary( string type);
    void setBottomBoundary(string type);
    void setLeftBoundary(string type);
    void setRightBoundary(string type);
    void setBoundaryNeighbours(int v);
    void setBoundaryNeighbours();
    void updateBoundaryVariables(int v);

    // Functions to create output file in VTK format
    void writeVTK(string fileName);

    void ResetMap_OutFlow(); // Function to Reset the Map to Outflow Boundaries.
    void updateOutFlowBoundary(int u, int v); // Function to update the map of Outflow Boundaries.
    void updateCellMarker(int v, int m);
    
    void setVanderMandMatrix(); // Function to set the VanderMand and its Inverse

    // Functions to handle Moment Limiter
    void computeMoments(int v, int m, int cm);
    void limitMoments(int m, int modifiedm, int cm, unsigned Index);
    void convertMomentToVariable(int m, int v, int cm);

    // Function to Reset Cell Centered Variables 
    void ResetVariables_CellCentered(int v, double value = 0.0);
    
    // Operators on the field.
    void delByDelX(int v, int vDash,int conserVar, string fluxType, int fluxVariable);

    void delByDelY(int v, int vDash, int conserVar, string fluxType, int fluxVariable);

    void delByDelX(int v, int vDash, string fluxType);

    void delByDelY(int v, int vDash, string fluxType);

    // Functions to apply linear operations on the variables.
    void axpy(double a, int x, int y);
    void scal(double a, int x);
    void setFunctionsForVariables(double a, int x, function<void(double, double*, unsigned, unsigned, double*)> f, int z); 
    void setFunctionsForVariables(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, unsigned, double*)> f, int z); 
    void setFunctionsForVariables(double a, int w, double b, int x, double c, int y, function<void(double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z);
    void setFunctionsForVariables(double t, int a, double u, int b, double v, int c, double x, int d, function<void(double, double*, double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z); 
     
    // Functions to apply operations only on variables stored at Boundary
    void setFunctionsForBoundaryVariables(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, unsigned, double*)> f, int z); 
    void setFunctionsForBoundaryVariables(double a, int w, double b, int x, double c, int y, function<void(double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z); 
    // Fucntions to compute cell centered values
    void setFunctionsForVariablesCellCentered(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, double*)> f, int z); 
    void setFunctionsForCellCenterVariablesfromDomainVariables(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, unsigned, double*, unsigned)> f, int z); 
    // Functions to give the information about the error.
    double l2Norm(int v1, int v2);

    // Function to check Positivity of data stored in cell .
    void checkPositivity(int v, int cm, string level);
    void resetPositivity();


    // Reworked Methods to update Boundary ,considering the enitre system of equation rather than individual variables

    void addConservativeVariables(int v);
    void addConservativeVariables(vector<int> V);

    void updateBoundary(double time);
    
    // Misc. Functions on for global field data 
    double FindMax(int v);
    double FindMindx();
    void FindMindx(int v);
    void FindTimestep(int dt, int dx, int U, double CFL);
    double FindMindt(int dt);
};

#endif
