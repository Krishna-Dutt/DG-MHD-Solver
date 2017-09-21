/**
 * @file DG_Field_2d.cpp
 * @Synopsis  This is the file for the class `Field`, which the data-structure for storing DG Fields
 * @author Kaushik Kulkarni
 * @version 1.0
 * @date 2017-02-18
 */

#include "../../includes/DG_Field_2d/DG_Field_2d.h"
#include "../../includes/DG_Element_2d/DG_Element_2d.h"
#include "../../includes/DG_BoundaryElement_2d/DG_BoundaryElement_2d.h"

#include "../../includes/Utilities/Inverse.h"
#include "../../includes/Utilities/Transpose.h"
#include "../../includes/Utilities/DerivativeMatrix.h"
#include "../../includes/Utilities/MassMatrix.h"
#include "../../includes/Utilities/FluxMatrix.h"
#include "../../includes/Utilities/VanderMandLegendre.h"




/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the constructor method which takes the following inputs. Once it gets the inputs, it starts
 * creating the Elements, which in turn starts allocating space to it, and also giving the inputs of domain information
 * to each and every element.
 *
 * @Param _nex The number of elements in the x-direction for the field
 * @Param _ney The number of elements in the y-direction for the field
 * @Param _N The order of interpolation of the polynomial which is to be used.
 * @Param _x1 The grid is a structured rectangular grid, and hence this corresponds to the x-coord of the bottom left
 * @Param _y1 This corresponds to the y-coordinate of the bottom left corner of the grid
 * corner of the grid.
 * @Param _x2 This corresponds to the x-coordinate of the top right corner of the grid
 * @Param _y2 This corresponds to the y- coordinate of the top right corner of the grid
 */
/* ----------------------------------------------------------------------------*/
DG_Field_2d::DG_Field_2d(int _nex, int _ney, int _N, double _x1, double _y1, double _x2, double _y2) {
    ne_x = _nex;
    ne_y = _ney;
    N = _N;
    x1 = _x1;
    x2 = _x2;
    y1 = _y1;
    y2 = _y2;
    
    /// Setting the grid, by setting the elements. The elements are set by providing their end points for the quads.
    elements.resize(ne_x);

    double x_curr,y_curr, dx = (x2-x1)/ne_x, dy = (y2-y1)/ne_y;
    x_curr = x1;
    y_curr = y1;

    // Setting up exponential grid (noo-uniform grid) along y direction
    double Beta_y, Beta_x, DeltaX, DeltaY, epsilon = 1e-10;
    DeltaY = (y2-y1);
    Beta_y = 1.1;

    DeltaX = (x2-x1);
    Beta_x = 1.2;
    dx = DeltaX * (Beta_x - 1.0 + epsilon)/(pow(Beta_x, ne_x) -1.0 + epsilon);

    for(int i=0; i<ne_x; i++){
        y_curr = y1;
        dy = DeltaY * (Beta_y - 1.0 + epsilon)/(pow(Beta_y, ne_y) -1.0 + epsilon);
        for(int j=0; j<ne_y; j++){
            if ( i == 0 || i == ne_x-1 || j == 0 || j == ne_y-1) {
                elements[i].push_back( new DG_BoundaryElement_2d(N, x_curr, y_curr, x_curr + dx, y_curr + dy) );
            }
            else
             {
                elements[i].push_back( new DG_Element_2d(N, x_curr, y_curr, x_curr + dx, y_curr + dy) );
            }

            y_curr += dy;
            dy = Beta_y*dy;
        }
        x_curr += dx;
        dx = Beta_x*dx;
    } // All the elements have been initialized.

    /// Setting the interaction between the elements by passing their neighboring elements addresses to each of the
    //elements.
    
    // Setting the top elements of each of the elements.
    for(int i = 0; i < ne_x; i++)
        for(int j=0; j < (ne_y-1); j++)
            elements[i][j]->setNeighboringElement('T', elements[i][j+1]);
        

    // Setting the right elements of each of the elements.
    for(int i = 0; i < (ne_x-1); i++)
        for(int j=0; j < (ne_y); j++)
            elements[i][j]->setNeighboringElement('R', elements[i+1][j]);

    // Setting the bottom elements of each of the elements.
    for(int i = 0; i < ne_x; i++)
        for(int j=1; j < ne_y; j++)
            elements[i][j]->setNeighboringElement('B', elements[i][j-1]);
    
    // Setting the left elements of each of the elements.
    for(int i = 1; i < ne_x; i++)
        for(int j=0; j < ne_y; j++)
            elements[i][j]->setNeighboringElement('L', elements[i-1][j]);
    // All the neighboring elements have been set, except for the elements at the boundary.
    
    /// Computing and passing the mass matrices, derivative and flux matrices to each and every element.
    /// The main use of this step is to ensure that each and every elements is not computing the same matrices.
    
    /// Assigning spaces to those matrices.
    double *massMatrix = new double[(N+1)*(N+1)*(N+1)*(N+1)];
    double *derivativeMatrix_x = new double[(N+1)*(N+1)*(N+1)*(N+1)] ; /// This is the matrix which is laid out in 1-D, this would help us to find t    he $\frac{d}{dx}$ of any term. 
    double *derivativeMatrix_y = new double[(N+1)*(N+1)*(N+1)*(N+1)] ; /// This is the matrix which is laid out in 1-D, this would help us to find t    he $\frac{d}{dy}$ of any term.
    double* fluxMatrix_top = new double[(N+1)*(N+1)*(N+1)*(N+1)] ; /// This is the flux matrix for the top edge.
    double* fluxMatrix_right = new double[(N+1)*(N+1)*(N+1)*(N+1)] ; // The Flux Matrix for the right edge.
    double* fluxMatrix_bottom = new double[(N+1)*(N+1)*(N+1)*(N+1)] ; /// This would be the flux term for the the bottom edge.
    double* fluxMatrix_left = new double[(N+1)*(N+1)*(N+1)*(N+1)] ; /// The Flux matrix for the left edge.
    double* massInverse = new double[(N+1)*(N+1)*(N+1)*(N+1)] ;
    double *TansposederivativeMatrix_x = new double[(N+1)*(N+1)*(N+1)*(N+1)] ;  
    double *TransposederivativeMatrix_y = new double[(N+1)*(N+1)*(N+1)*(N+1)] ;

    
    /// Calling functions to compute the entries of the matrix.
    
    twoDMassMatrix(massMatrix, N);
    inverse(massMatrix,massInverse,(N+1)*(N+1));
    twoDDerivativeMatrixX(derivativeMatrix_x, N);
    twoDDerivativeMatrixY(derivativeMatrix_y, N);
    transpose(derivativeMatrix_x, TansposederivativeMatrix_x, (N+1)*(N+1));
    transpose(derivativeMatrix_y, TransposederivativeMatrix_y, (N+1)*(N+1));
    
    twoDFluxMatrix3(fluxMatrix_top, N);
    twoDFluxMatrix2(fluxMatrix_right, N);
    twoDFluxMatrix4(fluxMatrix_left, N);
    twoDFluxMatrix1(fluxMatrix_bottom, N);

    /// Assigning the computed matrices to each and every element.
    for(int i = 0; i < ne_x; i++)
        for(int j=0; j < ne_y; j++){
            elements[i][j]->setMassMatrix(massMatrix);
            elements[i][j]->setInverseMassMatrix(massInverse);
            elements[i][j]->setderivateMatrix_x(derivativeMatrix_x);
            elements[i][j]->setderivateMatrix_y(derivativeMatrix_y);
            elements[i][j]->setTopFluxMatrix(fluxMatrix_top);
            elements[i][j]->setRightFluxMatrix(fluxMatrix_right);
            elements[i][j]->setLeftFluxMatrix(fluxMatrix_left);
            elements[i][j]->setBottomFluxMatrix(fluxMatrix_bottom);

            elements[i][j]->setTransposederivateMatrix_x(TansposederivativeMatrix_x);
            elements[i][j]->setTransposederivateMatrix_y(TransposederivativeMatrix_y);

        }

    
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function which finds and sets the VanderMandMatrix and its Inverse.
 * 
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setVanderMandMatrix() {
  double *vanderMand = new double [(N+1)*(N+1)*(N+1)*(N+1)];
  double *inverseVanderMand = new double [(N+1)*(N+1)*(N+1)*(N+1)];

  twoDVanderMandLegendre(vanderMand, N);
  inverse(vanderMand, inverseVanderMand, (N+1)*(N+1));

  for (int i=0; i < ne_x; ++i)
    for(int j=0; j < ne_y; ++j) {
     elements[i][j]->setVanderMandMatrix(vanderMand);
     elements[i][j]->setInverseVanderMandMatrix(inverseVanderMand);
    }

  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the destructor method  which deallocates and destroys the current object.
 * 
*/
/* ----------------------------------------------------------------------------*/
DG_Field_2d::~DG_Field_2d() {
  // To explicitly deallocate memory associated with the DG Matrices, shared by all elements
  elements[0][0]->Destroy_Matrices();
  for(int i = 0; i < ne_x; i++)
    for (int j = 0; j < ne_y; j++) {
      delete elements[i][j];
    }
}

void DG_Field_2d::setBoundaryConditions(string type) {
    if(type == "periodic") {
        // Setting the boundary for the top elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][ne_y-1]->setNeighboringElement('T', elements[i][0]);

        // Setting the boundary for the right elements.
        for(int j=0; j < (ne_y); j++)
            elements[ne_x-1][j]->setNeighboringElement('R', elements[0][j]);

        // Setting the boundary for the bottom elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][0]->setNeighboringElement('B', elements[i][ne_y-1]);

        // Setting the boundary for the left elements..
        for(int j=0; j < ne_y; j++)
            elements[0][j]->setNeighboringElement('L', elements[ne_x-1][j]);
    }
    else if(type == "periodicX") {
        
        // Setting the boundary for the right elements.
        for(int j=0; j < (ne_y); j++)
            elements[ne_x-1][j]->setNeighboringElement('R', elements[0][j]);

        // Setting the boundary for the left elements..
        for(int j=0; j < ne_y; j++)
            elements[0][j]->setNeighboringElement('L', elements[ne_x-1][j]);
    }
    else if(type == "periodicY") {
        // Setting the boundary for the top elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][ne_y-1]->setNeighboringElement('T', elements[i][0]);

        // Setting the boundary for the bottom elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][0]->setNeighboringElement('B', elements[i][ne_y-1]);
    }
    
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary type for a variable at Top Boundary.
 *
 * @Param v This is a string which defines the variable name.
 * @Param type This is a string which defines the boundary type.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setTopBoundary(string v, string type) {
    if(type == "periodic") {
        // Setting the boundary for the top elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][ne_y-1]->setNeighboringElement('T', elements[i][0]);


        // Setting the boundary for the bottom elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][0]->setNeighboringElement('B', elements[i][ne_y-1]);

    }
    else if(type == "neumann") {
        // Keep the neighbour cell pointer as this pointer
        
        
    }
    else if(type == "dirichlet") {
        // Keep the neighbour cell pointer as this pointer
       
    }

    else if(type == "reflective") {
        for(int i = 0; i < ne_x; i++)
            for(int j=0; j<=N; ++j) {
                *(elements[i][ne_y-1]->neighboringTop[v][j]) = -*(elements[i][ne_y-1]->boundaryTop[v][j]);
            }

    }

    for(int i = 0; i < ne_x; i++)
            elements[i][ne_y-1]->assignBoundary(v, type,'T');

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary type for a variable at Bottom Boundary.
 *
 * @Param v This is a string which defines the variable name.
 * @Param type This is a string which defines the boundary type.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setBottomBoundary(string v, string type) {
    if(type == "periodic") {
        // Setting the boundary for the top elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][ne_y-1]->setNeighboringElement('T', elements[i][0]);


        // Setting the boundary for the bottom elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][0]->setNeighboringElement('B', elements[i][ne_y-1]);

    }
    else if(type == "neumann") {
        // Keep the neighbour cell pointer as this pointer
        
        
    }
    else if(type == "dirichlet") {
        // Keep the neighbour cell pointer as this pointer
       
    }

    else if(type == "reflective") {
        for(int i = 0; i < ne_x; i++)
            for(int j=0; j<=N; ++j) {
                *(elements[i][0]->neighboringBottom[v][j]) = -*(elements[i][0]->boundaryBottom[v][j]);
            }

    }

    for(int i = 0; i < ne_x; i++)
            elements[i][0]->assignBoundary(v, type,'B');

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary type for a variable at Left Boundary.
 *
 * @Param v This is a string which defines the variable name.
 * @Param type This is a string which defines the boundary type.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setLeftBoundary(string v, string type) {
    if(type == "periodic") {
        // Setting the boundary for the right elements.
        for(int j=0; j < (ne_y); j++)
            elements[ne_x-1][j]->setNeighboringElement('R', elements[0][j]);

        // Setting the boundary for the left elements..
        for(int j=0; j < ne_y; j++)
            elements[0][j]->setNeighboringElement('L', elements[ne_x-1][j]);

    }
    else if(type == "neumann") {
        // Keep the neighbour cell pointer as this pointer
        
        
    }
    else if(type == "dirichlet") {
        // Keep the neighbour cell pointer as this pointer
       
    }

    else if(type == "reflective") {
        for(int i = 0; i < ne_y; i++)
            for(int j=0; j<=N; ++j) {
                *(elements[0][i]->neighboringLeft[v][j]) = -*(elements[0][i]->boundaryLeft[v][j]);
            }

    }

    for(int i = 0; i < ne_y; i++)
            elements[0][i]->assignBoundary(v, type,'L');

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary type for a variable at Right Boundary.
 *
 * @Param v This is a string which defines the variable name.
 * @Param type This is a string which defines the boundary type.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setRightBoundary(string v, string type) {
    if(type == "periodic") {
        // Setting the boundary for the right elements.
        for(int j=0; j < (ne_y); j++)
            elements[ne_x-1][j]->setNeighboringElement('R', elements[0][j]);

        // Setting the boundary for the left elements..
        for(int j=0; j < ne_y; j++)
            elements[0][j]->setNeighboringElement('L', elements[ne_x-1][j]);

    }
    else if(type == "neumann") {
        // Keep the neighbour cell pointer as this pointer
        
        
    }
    else if(type == "dirichlet") {
        // Keep the neighbour cell pointer as this pointer
       
    }

    else if(type == "reflective") {
        for(int i = 0; i < ne_y; i++)
            for(int j=0; j<=N; ++j) {
                *(elements[ne_x-1][i]->neighboringRight[v][j]) = -*(elements[ne_x-1][i]->boundaryRight[v][j]);
            }

    }

    for(int i = 0; i < ne_y; i++)
            elements[ne_x-1][i]->assignBoundary(v, type,'R');

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary  for a variable for Boundary Cells .
 *
 * @Param v This is a string which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setBoundaryNeighbours(string v) {
        // Setting the boundary for the top elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][ne_y-1]->setVariableNeighbors(v);


        // Setting the boundary for the bottom elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][0]->setVariableNeighbors(v);

        // Setting the boundary for the right elements.
        for(int j=0; j < (ne_y); j++)
            elements[ne_x-1][j]->setVariableNeighbors(v);

        // Setting the boundary for the left elements..
        for(int j=0; j < ne_y; j++)
            elements[0][j]->setVariableNeighbors(v);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates the given Variable at the Boundary .
 *
 * @Param v This is a string which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::updateBoundaryVariables(string v) {
        // Setting the boundary for the top elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][ne_y-1]->updateBoundaryVariables(v);


        // Setting the boundary for the bottom elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][0]->updateBoundaryVariables(v);

        // Setting the boundary for the right elements.
        for(int j=0; j < (ne_y); j++)
            elements[ne_x-1][j]->updateBoundaryVariables(v);

        // Setting the boundary for the left elements..
        for(int j=0; j < ne_y; j++)
            elements[0][j]->updateBoundaryVariables(v);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function adds the variable to each and every element. In this the boundary points are not
 * specifically stored in each element.
 *
 * @Param v This is a string which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::addVariable_withBounary(string v) {
    
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_withBoundary(v); // Adding the variable for the (i, j) th element.
       }
   }
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->setVariableNeighbors(v); // This is essential so that the addresses of the neighbors are stored in each and every element.
       }
   }
   variableNames.push_back(v);// What about vector variableWithBoundaryInfo, why is it there ??
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function adds the variable, required only at the boundary,to each and every element. In this the boundary points are not
 * specifically stored in each element.
 *
 * @Param v This is a string which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::addVariable_onlyBounary(string v) {
    
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_onlyBoundary(v); // Adding the variable for the (i, j) th element.
       }
   }
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->setVariableNeighbors(v); // This is essential so that the addresses of the neighbors are stored in each and every element.
       }
   }
  // variableNames.push_back(v); // Is this required, or rather, would it affect the solver in some way ; Shouldn't bee included in variablenames, as there are being written in output vtk file !!
   variableOnlyAtBoundary.push_back(v);
    return ;
}




/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is similar to DG_Field_2d::addVariable_withBoundary. Just the boundary points are not specificaly
 * stored for this variable.
 *
 * @Param v This is the name of the variable which is to be added.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::addVariable_withoutBounary(string v) {
    
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_withoutBoundary(v); // Adding the variable for the (i, j) th element.
       }
   }
  // variableNames.push_back(v); // These are support variables, Hence not required in output vtk files.
   return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is similar to DG_Field_2d::addVariable_withoutBoundary.Just only one value is store at for each cell,
 *  could be considered as cell centered value.
 *
 * @Param v This is the name of the variable which is to be added.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::addVariable_CellCentered(string v) {
    
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_CellCentered(v); // Adding the variable for the (i, j) th element.
       }
   }
   // variableNames.push_back(v);
   return ;
}


void DG_Field_2d::initializeVariable(string v, function<double(double, double)> f) {
    

   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
            elements[i][j]->initializeVariable(v, f); // Initializing the corresponding element by passing the same parameters to it.
       }
   }

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function in order to Reset the Cell Centered Variables to given value.
 *
 * @Param v The variable to be reset.
 * @Param value The value to which the variable is to be reset.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::ResetVariables_CellCentered(string v, double value) {
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
            elements[i][j]->ResetVariables_CellCentered(v, value);
       }
   }

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function in order to Reset the Map to OutFlow Boundaries.
 *
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::ResetMap_OutFlow() {
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
            elements[i][j]->ResetMap_OutFlow();
       }
   }

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to update the Map to OutFlow Boundaries.
 * 
 * @Param u This is the velocity in x direction.
 * @Param v This is the velocity in y direction.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::updateOutFlowBoundary(string u, string v) {
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
            elements[i][j]->updateOutFlowBoundary(u, v);
       }
   }

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to update the Cell Markers for Shock Detection.
 * 
 * @Param v This is the quantity used for detecting Shocks/Discontinuities.
 * @Param m This is the variable used to store the value of cell marker.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::updateCellMarker(string v, string m) {
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
            elements[i][j]->updateCellMarker(v, m);
       }
   }

    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to compute the Moments for given variable.
 * 
 * @Param v This is the quantity whose moments are to be calculated.
 * @Param m This is the variable used to store the moments.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::computeMoments(string v, string m) {
  for(int i=0; i < ne_x; ++i)
    for(int j=0; j < ne_y; ++j) {
      elements[i][j]->computeMoments(v, m);
    }

  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to limit the Moments for given variable.
 * 
 * @Param m This is gives the moments to be limited.
 * @Param modifiedm This is the variable to store the modified moments.
 * @Param cm This the cell marker used to identified troubled cells.
 * @Param Index This is the index at which to start limiting .
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::limitMoments(string m, string modifiedm, string cm, unsigned Index) {
  for(int i=0; i < ne_x; ++i)
    for(int j=0; j < ne_y; ++j) {
      elements[i][j]->limitMoments(m, modifiedm, cm, Index);
    }

  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to obtain value of the variable from the moments.
 * 
 * @Param m This is gives the moments of the variable.
 * @Param v This is the variable /quantity to be obtained.
 * @Param cm This the cell marker used to identified troubled cells.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::convertMomentToVariable(string m, string v, string cm) {
  for(int i=0; i < ne_x; ++i)
    for(int j=0; j < ne_y; ++j) {
      elements[i][j]->convertMomentToVariable(m, v, cm);
    }

  return ;
}



/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function in order to write the data in the form of VTK file.
 *
 * @Param fileName This is the string fileName with which the file is to be saved.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::writeVTK(string fileName){
    ofstream pFile;
    pFile.open(fileName);
    
    int i, j, k, k1, k2;

    // Printing the preamble for the .vtk file.
    pFile << "# vtk DataFile Version 3.0\nNavier Stokes DG\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    // The information of the number of points.
    pFile << "POINTS\t" << (N+1)*(N+1)*ne_x*ne_y << "\tdouble\n";

    // Writing the point co-ordinates.
    for ( j = 0; j < ne_y; j++ )
        for ( i = 0; i < ne_x; i++ )
            for( k = 0; k < (N+1)*(N+1); k++ )
                pFile << elements[i][j]->X[k] << "\t" << elements[i][j]->Y[k] <<"\t"<<0 <<endl;

    pFile << "\n\n";

    // Specifying the information about the CELLS.
    pFile << "CELLS\t" << (N*N*ne_x*ne_y) <<"\t" << 5*(N*N*ne_x*ne_y) << endl;

    // Writing information about the structure of the cells.
    for ( i = 0; i < ne_y; i++ ) {
        for ( j = 0; j < ne_x; j++ ) {
            for( k1 = 0; k1 < N; k1++ ) {
                for ( k2 = 0; k2 < N; k2++ ) {
                    k   =   (i*ne_x+j)*(N+1)*(N+1) +   k1*(N+1)    +   k2;
                    pFile << 4 << "\t" << k << "\t" << k+1 << "\t" << k+N+2 << "\t" << k+N+1 << endl;
                }
            }
        }
    }
    pFile << "\n\n";

    // Specifying the information about the CELL TYPES.
    pFile << "CELL_TYPES " << (N*N*ne_x*ne_y) << endl;

    // `9` is the CELL TYPE CODE for specifying that it is a quad.
    for ( i = 0; i < (N*N*ne_x*ne_y); i++)
        pFile << "9\n";
    pFile << "\n\n";

    // Specifying the information about the values of the scalars.
    
    pFile << "POINT_DATA\t"<< (N+1)*(N+1)*ne_x*ne_y <<"\n";
    
    int noOfVars = variableNames.size(); // Getting the number of variables
    double* currentVariable;

    for(k1=0; k1<noOfVars; k1++) {
        pFile << "SCALARS\t"<< variableNames[k1] <<"\tdouble\nLOOKUP_TABLE default\n";
        
        // Writing the value of the POINT_DATA, for the variable[variableNames[k1]] 
        for ( j = 0; j < ne_y; j++ ){
            for ( i = 0; i < ne_x; i++ ) {
                currentVariable = elements[i][j]->variable[variableNames[k1]];
                for( k = 0; k < (N+1)*(N+1); k++ ) {
                    pFile << currentVariable[k] << endl;
                }
            }
        }
    }
    pFile.close(); // Closing the file.
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to operate the partial derivative of the variable w.r.t. x.
 *
 * @Param v The variable which is to be differentiated
 * @Param vDash The variable in which the differentiated value is to be stored.
 * @Param conserVar The correpsonding conservative/dependent variable.
 * @Param fluxType The numerical flux type which is to be implemented while computing the derivative.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::delByDelX(string v, string vDash, string conserVar, string fluxType, string fluxVariable = "") {
    for(int i = 0; i < ne_x; i++ )
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->delByDelX(v, vDash, conserVar, fluxType, fluxVariable);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to operate the partial derivative of the variable w.r.t. y.
 *
 * @Param v The variable which is to be differentiated
 * @Param vDash The variable in which the differentiated value is to be stored.
 * @Param conserVar The corresponding conservative/Dependent Variable.
 * @Param fluxType The numerical flux type which is to be implemented while computing the derivative.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::delByDelY(string v, string vDash, string conserVar, string fluxType, string fluxVariable = "") {
    for(int i = 0; i < ne_x; i++ )
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->delByDelY(v, vDash, conserVar, fluxType, fluxVariable);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to operate the partial derivative of the variable w.r.t. x.
 *
 * @Param v The variable which is to be differentiated
 * @Param vDash The variable in which the differentiated value is to be stored.
 * @Param fluxType The numerical flux type which is to be implemented while computing the derivative.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::delByDelX(string v, string vDash, string fluxType) {
    for(int i = 0; i < ne_x; i++ )
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->delByDelX(v, vDash, fluxType);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to operate the partial derivative of the variable w.r.t. y.
 *
 * @Param v The variable which is to be differentiated
 * @Param vDash The variable in which the differentiated value is to be stored.
 * @Param fluxType The numerical flux type which is to be implemented while computing the derivative.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::delByDelY(string v, string vDash, string fluxType) {
    for(int i = 0; i < ne_x; i++ )
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->delByDelY(v, vDash, fluxType);

    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @brief      This is used to extern the function `saxpy` to the class DG_Field_2d. $\mathbf{y} \mapsto a\mathbf{x} + \mathbf{y}$
 *
 * @param[in]  a     The coefficient `a` of `x`
 * @param[in]  x     The column vector `x`
 * @param[in]  y     The column vector `y`
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::axpy(double a, string x, string y) {

    for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->axpy(a, x, y);
    
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @brief      This is used to extern the function `sscal` to the class DG_Field_2d. $\mathbf{x} \mapsto a\mathbf{x}$
 *
 * @param[in]  a     The coefficient `a` of `x`
 * @param[in]  x     The column vector `x`
 * @param[in]  y     The column vector `y`
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::scal(double a, string x) {
    for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->scal(a, x);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable z to f(x, y).
 *
 * @Param x The first parameter of the function.
 * @Param y The second parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForVariables(string x, string y, function<double(double, double)> f, string z) {
    for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->setFunctionsForVariables(x, y, f, z);
    return;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable z to f(w, x, y).
 *
 * @Param w The first parameter of the function
 * @Param x The second parameter of the function.
 * @Param y The third parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForVariables(string w, string x, string y, function<double(double, double, double)> f, string z) {
    for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->setFunctionsForVariables(w, x, y, f, z);
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable z to f(a, b, c, d).
 *
 * @Param a The first parameter of the function
 * @Param b The second parameter of the function.
 * @Param c The third parameter of the function.
 * @Param d The third parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForVariables(string a, string b, string c, string d, function<double(double, double, double, double)> f, string z) {
    for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->setFunctionsForVariables(a, b, c, d, f, z);
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of Boundary variable z to f(x, y).
 *
 * @Param x The first parameter of the function.
 * @Param y The second parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForBoundaryVariables(string x, string y, function<double(double, double)> f, string z) {
    for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->setFunctionsForBoundaryVariables(x, y, f, z);
    return;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of Boundary variable z to f(w, x, y).
 *
 * @Param w The first parameter of the function
 * @Param x The second parameter of the function.
 * @Param y The third parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForBoundaryVariables(string w, string x, string y, function<double(double, double, double)> f, string z) {
    for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->setFunctionsForBoundaryVariables(w, x, y, f, z);
    return ;
}



double DG_Field_2d::l2Norm(string v1, string v2) {
    double norm = 0.0;
    double elementNorm;
    for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++){
            elementNorm = elements[i][j]->l2Norm(v1, v2);
            norm += (elementNorm*elementNorm);
        }

    norm = sqrt(norm/(ne_x*ne_y));    
    return norm;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function checks the positivity of data stored in cell.
 *
 * @Param v  This is the name of the variable whose positivity is to be checked.
 * @Param cm This is the cell marker used to identify troubled cells.
 * @Param level This string is used to identify the level of limiting rquired.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::checkPositivity(string v, string cm, string level) {
     for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->checkPositivity(v, cm, level);
    return;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function resets the positivity markers in cell.
 *
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::resetPositivity() {
     for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->resetPositivity();
    return;
}

