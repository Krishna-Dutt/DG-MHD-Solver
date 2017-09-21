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
    vector<string> variableNames; // This is stores all the variables which have been added to the field.
    vector<string> variablesWithBoundaryInfo; // This stores all the variables whose boundary info. is also known.

    vector<string> variableOnlyAtBoundary; // This stores all the variables who are required only at the Boundaries.

    DG_Field_2d(int _nex, int _ney, int _N, double _x1, double _y1, double _x2, double _y2);
    ~DG_Field_2d();

    void setSystem(string S);

    void addVariable_withBounary(string v);
    void addVariable_withoutBounary(string v);
    void addVariable_onlyBounary( string v);
    void addVariable_CellCentered(string v);
    void initializeVariable(string v, function<double(double, double)>);

    // Functions to set Boundary Conditions
    void setBoundaryConditions(string type);
    void setBoundaryConditions(string v, string type, string b);
    void setTopBoundary(string v, string type);
    void setBottomBoundary(string v, string type);
    void setLeftBoundary(string v, string type);
    void setRightBoundary(string v, string type);
    void setBoundaryNeighbours(string v);
    void updateBoundaryVariables(string v);

    // Functions to create output file in VTK format
    void writeVTK(string fileName);

    void ResetMap_OutFlow(); // Function to Reset the Map to Outflow Boundaries.
    void updateOutFlowBoundary(string u, string v); // Function to update the map of Outflow Boundaries.
    void updateCellMarker(string v, string m);
    
    void setVanderMandMatrix(); // Function to set the VanderMand and its Inverse

    // Functions to handle Moment Limiter
    void computeMoments(string v, string m);
    void limitMoments(string m, string modifiedm, string cm, unsigned Index);
    void convertMomentToVariable(string m, string v, string cm);

    // Function to Reset Cell Centered Variables 
    void ResetVariables_CellCentered(string v, double value = 0.0);
    
    // Operators on the field.
    void delByDelX(string v, string vDash,string conserVar, string fluxType, string fluxVariable);

    void delByDelY(string v, string vDash, string conserVar, string fluxType, string fluxVariable);

    void delByDelX(string v, string vDash, string fluxType);

    void delByDelY(string v, string vDash, string fluxType);

    // Functions to apply linear operations on the variables.
    void axpy(double a, string x, string y);
    void scal(double a, string x);
    void setFunctionsForVariables(string x, string y, function<double(double, double)>, string z); 
    void setFunctionsForVariables(string w, string x, string y, function<double(double, double, double)>, string z); 
    // Functions to apply operations only on variables stored at Boundary
    void setFunctionsForBoundaryVariables(string x, string y, function<double(double, double)>, string z); 
    void setFunctionsForBoundaryVariables(string w, string x, string y, function<double(double, double, double)>, string z); 
    void setFunctionsForVariables(string a, string b, string c, string d, function<double(double, double, double, double)>, string z); 

    // Functions to give the information about the error.
    double l2Norm(string v1, string v2);

    // Function to check Positivity of data stored in cell .
    void checkPositivity(string v, string cm, string level);
    void resetPositivity();


    // Reworked Methods to update Boundary ,considering the enitre system of equation rather than individual variables

    void addConservativeVariables(string v);
    void addConservativeVariables(vector<string> V);

    void updateBoundary();
};

#endif
