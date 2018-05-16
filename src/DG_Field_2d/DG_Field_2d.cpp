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
    system = "EULER";
    
    /// Setting the grid, by setting the elements. The elements are set by providing their end points for the quads.
    elements.resize(ne_x);

    double x_curr,y_curr, dx2, dx = (x2-x1)/ne_x, dy = (y2-y1)/ne_y;
    x_curr = x1;
    y_curr = y1;

    // Setting up exponential grid (noo-uniform grid) along y direction
    double Beta_y, Beta_x1, Beta_x2, DeltaX1, DeltaX2, DeltaY, epsilon = 1e-10;
    // Setting up Hyperbolic grids along x-direction about scale_b X (x2-x1)
    scale_b = 1;//0.5;
    DeltaY = (y2-y1);
    Beta_y = 1.1;
    
    DeltaX1 = (x2-x1)*scale_b;
    Beta_x1 = 1.005;//1/1.00004;
    dx = DeltaX1 * (Beta_x1 - 1.0 + epsilon)/(pow(Beta_x1, ne_x*scale_b) -1.0 + epsilon);
    
    for(int i=0; i<ne_x*scale_b; i++){
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
        dx = Beta_x1*dx;
    } // All the elements have been initialized.

    DeltaX2 = (x2-x1)*(1-scale_b);
    Beta_x2 = 1.00004;
    //dx = DeltaX2 * (Beta_x2 - 1.0 + epsilon)/(pow(Beta_x2, ne_x*(1-scale_b)) -1.0 + epsilon);

    for(int i=ne_x*scale_b; i<ne_x; i++){
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
        //dx = Beta_x2*dx;
    } // All the elements have been initialized.

    /// Setting the interaction between the elements by passing their neighboring elements addresses to each of the
    //elements.
    
    // Setting the top elements of each of the elements.
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j=0; j < (ne_y-1); j++)
            elements[i][j]->setNeighboringElement('T', elements[i][j+1]);
        

    // Setting the right elements of each of the elements.
    #pragma omp parallel for private(j)
    for(int i = 0; i < (ne_x-1); i++)
        for(j=0; j < (ne_y); j++)
            elements[i][j]->setNeighboringElement('R', elements[i+1][j]);

    // Setting the bottom elements of each of the elements.
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j=1; j < ne_y; j++)
            elements[i][j]->setNeighboringElement('B', elements[i][j-1]);
    
    // Setting the left elements of each of the elements.
    #pragma omp parallel for private(j)
    for(int i = 1; i < ne_x; i++)
        for(j=0; j < ne_y; j++)
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
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j=0; j < ne_y; j++){
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
    
    variableNames.resize(0);
    variablesWithBoundaryInfo.resize(0);
    variableOnlyAtBoundary.resize(0);
    domainVariable.resize(0);
    cellcenterVariable.resize(0);
    
    PositivityMarker = new bool[ne_x*ne_y];
    for(int i =0 ; i < ne_x*ne_y ; ++i) 
       PositivityMarker[i] = true;
    
    RightEigenMatrix = NULL;
    LeftEigenMatrix = NULL;
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

  // Printing both matrices 
 /* cout << "Vander Mand Matrix  :: \n";
  for ( int i=0; i < (N+1)*(N+1) ; ++i) {
      for (j=0; j < (N+1)*(N+1); ++j) {
          cout << vanderMand[i*(N+1)*(N+1) + j] << " : ";
      }
      cout << "\n";
  }

  cout << "Inverse Vander Mand Matrix  :: \n";
  for ( int i=0; i < (N+1)*(N+1) ; ++i) {
      for (j=0; j < (N+1)*(N+1); ++j) {
          cout << inverseVanderMand[i*(N+1)*(N+1) + j] << " : ";
      }
      cout << "\n";
  }*/
  int j;
  #pragma omp parallel for private(j)
  for (int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j) {
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
  delete[] PositivityMarker;

  if(RightEigenMatrix != NULL) {
      delete[] RightEigenMatrix;
      delete[] LeftEigenMatrix;
  }

  for(int i=0; i < domainVariable.size(); ++i) {
      delete[] domainVariable[i];
  }
  for(int i=0; i < cellcenterVariable.size(); ++i) {
      delete[] cellcenterVariable[i];
  }
  for(int i = 0; i < ne_x; i++)
    for (int j = 0; j < ne_y; j++) {
      delete elements[i][j];
    }
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the System of Governing Euqaitons ot be solved..
 *
 * @Param S The governing equations
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setSystem(string S) {
    system = S;
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
    for (j = 0; j < ne_y; j++) {
      elements[i][j]->setSystem(S);
    }

    return ;
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
 * @Param type This is a string which defines the boundary type.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setTopBoundary(string type) {
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

   

    for(int i = 0; i < ne_x; i++)
            elements[i][ne_y-1]->assignBoundary(type,'T');

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary type for a variable at Bottom Boundary.
 *
 * @Param type This is a string which defines the boundary type.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setBottomBoundary(string type) {
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

   

    for(int i = 0; i < ne_x; i++)
            elements[i][0]->assignBoundary( type,'B');

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary type for a variable at Left Boundary.
 *
 * @Param type This is a string which defines the boundary type.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setLeftBoundary( string type) {
     //cout << "Calling Left Boundary !!" << endl;
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

    

    for(int i = 0; i < ne_y; i++)
            elements[0][i]->assignBoundary(type,'L');

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary type for a variable at Right Boundary.
 *
 * @Param type This is a string which defines the boundary type.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setRightBoundary( string type) {
   // cout << "Calling Right Boundary !!" << endl;
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


    for(int i = 0; i < ne_y; i++)
            elements[ne_x-1][i]->assignBoundary(type,'R');

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary  for a variable for Boundary Cells .
 *
 * @Param v This is a int which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setBoundaryNeighbours(int v) {
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
 * @Synopsis  This function sets the boundary  for all boundary variable for Boundary Cells .
 *
 * @Param v This is a int which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setBoundaryNeighbours() {
    for(int k=0; k <variablesWithBoundaryInfo.size(); ++k) {
        // Setting the boundary for the top elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][ne_y-1]->setVariableNeighbors(variablesWithBoundaryInfo[k]);


        // Setting the boundary for the bottom elements.
        for(int i = 0; i < ne_x; i++)
            elements[i][0]->setVariableNeighbors(variablesWithBoundaryInfo[k]);

        // Setting the boundary for the right elements.
        for(int j=0; j < (ne_y); j++)
            elements[ne_x-1][j]->setVariableNeighbors(variablesWithBoundaryInfo[k]);

        // Setting the boundary for the left elements..
        for(int j=0; j < ne_y; j++)
            elements[0][j]->setVariableNeighbors(variablesWithBoundaryInfo[k]);
    }

    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates the given Variable at the Boundary .
 *
 * @Param v This is a int which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::updateBoundaryVariables(int v) {
        // Setting the boundary for the top elements.
        #pragma omp parallel for 
        for(int i = 0; i < ne_x; i++)
            elements[i][ne_y-1]->updateBoundaryVariables(v);


        // Setting the boundary for the bottom elements.
        #pragma omp parallel for 
        for(int i = 0; i < ne_x; i++)
            elements[i][0]->updateBoundaryVariables(v);

        // Setting the boundary for the right elements.
        #pragma omp parallel for 
        for(int j=0; j < (ne_y); j++)
            elements[ne_x-1][j]->updateBoundaryVariables(v);

        // Setting the boundary for the left elements.
        #pragma omp parallel for 
        for(int j=0; j < ne_y; j++)
            elements[0][j]->updateBoundaryVariables(v);

    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function adds the variable, required only at the boundary,to each and every element. In this the boundary points are not
 * specifically stored in each element.
 *
 * @Param v This is a int which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
/*void DG_Field_2d::addVariable_onlyBounary(int v) {
    
   for (int i=0; i < ne_x; i++ ){
       for (j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_onlyBoundary(v); // Adding the variable for the (i, j) th element.
       }
   }
   for (int i=0; i < ne_x; i++ ){
       for (j=0; j<ne_y; j++) {
           elements[i][j]->setVariableNeighbors(v); // This is essential so that the addresses of the neighbors are stored in each and every element.
       }
   }
  // variableNames.push_back(v); // Is this required, or rather, would it affect the solver in some way ; Shouldn't bee included in variablenames, as there are being written in output vtk file !!
   variableOnlyAtBoundary.push_back(v);
    return ;
}*/


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function adds the variable to each and every element. In this the boundary points are not
 * specifically stored in each element.
 * 
 * @Param V The name of the variable begin added
 */
/* ----------------------------------------------------------------------------*/
int DG_Field_2d::addVariable_withBounary(string V) {
    double *newVariable = new double[(N+1)*(N+1)*ne_x*ne_y]; /// Allocating the space for the new variable which is to be created.
    int v;

   domainVariable.push_back(newVariable); /// Now assigning the same to the map. 
   v = domainVariable.size()-1;
    
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_withBoundary(v, newVariable + (N+1)*(N+1)*(ne_y*i + j) ); // Adding the variable for the (i, j) th element.
       }
   }
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->setVariableNeighbors(v); // This is essential so that the addresses of the neighbors are stored in each and every element.
       }
   }
   boundaryVariableNames[v] = V;
   variablesWithBoundaryInfo.push_back(v);
   variableNames.push_back(v);

   //scal(0.0, v);

    return v;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is similar to DG_Field_2d::addVariable_withBoundary. Just the boundary points are not specificaly
 * stored for this variable.
 *
 */
/* ----------------------------------------------------------------------------*/
int DG_Field_2d::addVariable_withoutBounary() {
   double *newVariable = new double[(N+1)*(N+1)*ne_x*ne_y]; /// Allocating the space for the new variable which is to be created.
   int v;
   
   domainVariable.push_back(newVariable); /// Now assigning the same to the map. 
   v = domainVariable.size()-1;

   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_withoutBoundary(v, newVariable + (N+1)*(N+1)*(ne_y*i + j) ); // Adding the variable for the (i, j) th element.
       }
   }
   
   variableNames.push_back(v);
   //scal(0.0, v);

   return v;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is similar to DG_Field_2d::addVariable_withoutBoundary.Just only one value is store at for each cell,
 *  could be considered as cell centered value.
 *
 */
/* ----------------------------------------------------------------------------*/
int DG_Field_2d::addVariable_CellCentered() {
   double *newVariable = new double[ne_x*ne_y]; /// Allocating the space for the new variable which is to be created.
   int v;

   cellcenterVariable.push_back(newVariable); /// Now assigning the same to the map. 
   v = cellcenterVariable.size()-1;

   return v;
}


void DG_Field_2d::initializeVariable(int v, function<double(double, double)> f) {
   int j; 
   #pragma omp parallel for private(j)
   for (int i=0; i < ne_x; i++ ){
       for (j=0; j<ne_y; j++) {
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
void DG_Field_2d::ResetVariables_CellCentered(int v, double value) {
   int j;
   #pragma omp parallel for private(j) 
   for (int i=0; i < ne_x; i++ ){
       for (j=0; j<ne_y; j++) {
            cellcenterVariable[v][i*ne_y + j] = value;
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
   int j;
   #pragma omp parallel for private(j) 
   for (int i=0; i < ne_x; i++ ){
       for (j=0; j<ne_y; j++) {
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
void DG_Field_2d::updateOutFlowBoundary(int u, int v) {
   int j;
   #pragma omp parallel for private(j)
   for (int i=0; i < ne_x; i++ ){
       for (j=0; j<ne_y; j++) {
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
void DG_Field_2d::updateCellMarker(int v, int m) {
   int j;
   #pragma omp parallel for private(j) 
   for (int i=0; i < ne_x; i++ ){
       for (j=0; j<ne_y; j++) {
            cellcenterVariable[m][i*ne_y + j] = elements[i][j]->updateCellMarker(v);
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
 * @Param cm This the cell marker used to identified troubled cells.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::computeMoments(int v, int m, int cm) {
  int j;
  #pragma omp parallel for private(j)    
  for(int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j) {
        if ( cellcenterVariable[cm][i*ne_y + j] && PositivityMarker[i*ne_y + j]) {
             elements[i][j]->computeMoments(v, m);
        }
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
void DG_Field_2d::limitMoments(int m, int modifiedm, int cm, unsigned Index) {
  int j;
  #pragma omp parallel for private(j)   
  for(int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j) {
        if ( cellcenterVariable[cm][i*ne_y + j] && PositivityMarker[i*ne_y + j]) {
            elements[i][j]->limitMoments(m, modifiedm, Index);
        }
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
void DG_Field_2d::convertMomentToVariable(int m, int v, int cm) {
  int j;
  #pragma omp parallel for private(j)  
  for(int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j) {
        if ( cellcenterVariable[cm][i*ne_y + j] && PositivityMarker[i*ne_y + j] ) {
            elements[i][j]->convertMomentToVariable(m, v);
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
 * @Param cm This the cell marker used to identified troubled cells.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::computeMoments(int *v, int *m, int cm, unsigned size) {
  int j;
  #pragma omp parallel for private(j)  
  for(int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j) {
        if ( cellcenterVariable[cm][i*ne_y + j] && PositivityMarker[i*ne_y + j]) {
             elements[i][j]->computeMoments(v, m, size);
        }
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
void DG_Field_2d::limitMoments(int *m, int *modifiedm, int cm, unsigned Index, unsigned size) {
  int j;
  #pragma omp parallel for private(j)  
  for(int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j) {
        if ( cellcenterVariable[cm][i*ne_y + j] && PositivityMarker[i*ne_y + j]) {
            elements[i][j]->limitMoments(m, modifiedm, Index, size);
        }
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
void DG_Field_2d::convertMomentToVariable(int *m, int *v, int cm, unsigned size) {
  int j;
  #pragma omp parallel for private(j)  
  for(int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j) {
        if ( cellcenterVariable[cm][i*ne_y + j] && PositivityMarker[i*ne_y + j] ) {
            elements[i][j]->convertMomentToVariable(m, v, size);
        }
    }

  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to limit the Moments for given variable, using Characteristic Limiter.
 * 
 * @Param array V array of indices of moments of Conservative Variables.
 * @Param C moment of characteristic variables, to store modified moments.
 * @Param cm This the cell marker used to identified troubled cells.
 * @Param Index Index correspoding to Characteristic variable.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::limitMoments(int *V, int *C, int cm, unsigned Index) {
  int j;
  #pragma omp parallel for private(j)  
  for(int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j)  {
        if ( cellcenterVariable[cm][i*ne_y + j] ) {
      elements[i][j]->limitMoments(V, C, Index);
      }
    }

  return ;
}



/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function which sets the Right and Left Eigen Matrix corresponding to directionalong flow
 * 
 * @Param dimension Dimension os system of Equations
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setEigenMatrices(unsigned _dimension) {
  
  Dimension = _dimension;
  RightEigenMatrix = new double[Dimension*Dimension*ne_x*ne_y];
  LeftEigenMatrix = new double[Dimension*Dimension*ne_x*ne_y];

  for (int i=0; i < ne_x; ++i)
    for(int j=0; j < ne_y; ++j) {
     elements[i][j]->setEigenMatrices(_dimension, (RightEigenMatrix + Dimension*Dimension*(i*ne_y + j)), (LeftEigenMatrix + Dimension*Dimension*(i*ne_y + j)));
    }

  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to find  the Right and Left Eigen Matrix corresponding to direction along flow
 * 
 * @Param array V array of strings corresponding to field variables required to compute eigne matrices.
 * @Param cm Cell Marker ,to identify troubled cells.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::findEigenMatrices(int *V, int cm) {
  int j;
  if ( system == "EULER") {
      #pragma omp parallel for private(j)
      for (int i=0; i < ne_x; ++i)
      for(j=0; j < ne_y; ++j)  {
        if ( cellcenterVariable[cm][i*ne_y + j] ) {
            elements[i][j]->findEigenMatricesEuler( V);
        }
      }

  }
  else if ( system == "MHD") {
      #pragma omp parallel for private(j)
      for (int i=0; i < ne_x; ++i)
      for(j=0; j < ne_y; ++j)  {
        if ( cellcenterVariable[cm][i*ne_y + j] ) {
            elements[i][j]->findEigenMatricesMHD(V);
        }
      }

  }
  
  return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to compute Characteristic Variables.
 * 
 * @Param array V Set of conservative variables..
 * @Param array C The Characteristic Variable.
 * @Param I Identifier for Characteristic Variable.
 * @Param cm Cell Marker to identify troubled cell.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::convertVariabletoCharacteristic(int *V, int *C, unsigned I, int cm) {
  int j;
  #pragma omp parallel for private(j)  
  for(int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j)  {
        if ( cellcenterVariable[cm][i*ne_y + j] ) {
              elements[i][j]->convertVariabletoCharacteristic(V, C, I);
        }
    }

  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to compute Conservative  Variables from Characteristic Variables.
 * 
 * @Param array C Set of characteristic variables..
 * @Param array V The Conservative Variable.
 * @Param I Identifier for conservative Variable.
 * @Param cm Cell Marker to identify troubled cell.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::convertCharacteristictoVariable(int *C, int *V, unsigned I, int cm) {
  int j;
  #pragma omp parallel for private(j)  
  for(int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j)  {
        if ( cellcenterVariable[cm][i*ne_y + j] ) {
             elements[i][j]->convertCharacteristictoVariable(C, V, I);
        }
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
    
    //int noOfVars = variableNames.size(); // Getting the number of variables
    double* currentVariable;

    for(map<int,string>::iterator itr = boundaryVariableNames.begin(); itr!= boundaryVariableNames.end(); ++itr) {
        pFile << "SCALARS\t"<< itr->second <<"\tdouble\nLOOKUP_TABLE default\n";
        
        // Writing the value of the POINT_DATA, for the variable
        for ( j = 0; j < ne_y; j++ ){
            for ( i = 0; i < ne_x; i++ ) {
                currentVariable = elements[i][j]->variable[itr->first];
                for( k = 0; k < (N+1)*(N+1); k++ ) {
                    pFile << currentVariable[k] << endl;
                }
            }
        }
    }

    /*for(k1=0; k1<domainVariable.size(); ++k1) {
        pFile << "SCALARS\t"<< domainVariable[k1] <<"\tdouble\nLOOKUP_TABLE default\n";
        for(i=0 ; i<(N+1)*(N+1)*ne_x*ne_y; ++i) {
            pFile << domainVariable[k1] << endl;
        }
    }*/

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
void DG_Field_2d::delByDelX(int v, int vDash, int conserVar, string fluxType, int fluxVariable = 999) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++ )
        for(j = 0; j < ne_y; j++)
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
void DG_Field_2d::delByDelY(int v, int vDash, int conserVar, string fluxType, int fluxVariable = 999) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++ )
        for(j = 0; j < ne_y; j++)
            elements[i][j]->delByDelY(v, vDash, conserVar, fluxType, fluxVariable);

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
void DG_Field_2d::delByDelX(int *V, int *VDash, int *ConserVar, string fluxType, int *fluxVariable, unsigned size) {
    int j; 
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++ )
        for(j = 0; j < ne_y; j++)
            elements[i][j]->delByDelX(V, VDash, ConserVar, fluxType, fluxVariable, size);

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
void DG_Field_2d::delByDelY(int *V, int *VDash, int *ConserVar, string fluxType, int *fluxVariable, unsigned size) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++ )
        for(j = 0; j < ne_y; j++)
            elements[i][j]->delByDelY(V, VDash, ConserVar, fluxType, fluxVariable, size);

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
void DG_Field_2d::delByDelX(int v, int vDash, string fluxType) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++ )
        for(j = 0; j < ne_y; j++)
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
void DG_Field_2d::delByDelY(int v, int vDash, string fluxType) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++ )
        for(j = 0; j < ne_y; j++)
            elements[i][j]->delByDelY(v, vDash, fluxType);

    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to compute the Divergence of Magnetic field.
 * 
 * @Param Bx This is Magnetic field in x direction.
 * @Param By This is Magnetic field in y direction.
 * @Param DeldotB This is used to store the computed locally constant Divergence of B.
*/
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::updateDivergenceB(int Bx, int By, int DeldotB) {
  int j;
  #pragma omp parallel for private(j)
  for(int i=0; i < ne_x; ++i)
    for(j=0; j < ne_y; ++j) {
        {
             elements[i][j]->updateDivergenceB(Bx, By, DeldotB);
        }
    }

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
void DG_Field_2d::axpy(double a, int x, int y) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            elements[i][j]->axpy(a, x, y);
    
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @brief      This is used to extern the function `sscal` to the class DG_Field_2d. $\mathbf{x} \mapsto a\mathbf{x}$
 *
 * @param[in]  a     The coefficient `a` of `x`
 * @param[in]  x     The column vector `x`
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::scal(double a, int x) {
    /*for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++) 
            elements[i][j]->scal(a, x);*/
    double *domVar = domainVariable[x];
    int j;
    if (!a) {
        #pragma omp parallel for 
        for(int i=0; i< (N+1)*(N+1)*ne_x*ne_y; ++i){
            domVar[i] = 0;
        }
        return ;
    }
    
    #pragma omp parallel for 
    for(int i=0 ; i< (N+1)*(N+1)*ne_x*ne_y ; ++i) {
        domVar[i] = a * domVar[i];
    }

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis     This is used to set a constasnt value to given vector
 *
 * @param a     The constant a
 * @param v     The column vector v
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setConstant(double a, int v) {
    /*for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++) 
            elements[i][j]->scal(a, x);*/
    double *domVar = domainVariable[v];
    
    #pragma omp parallel for 
    for(int i=0 ; i< (N+1)*(N+1)*ne_x*ne_y ; ++i) {
        domVar[i] = a;
    }

    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable z to f(x).
 *
 * @Param a Scaling value for x
 * @Param x The first parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForVariables(double a, int x, function<void(double, double*, unsigned, unsigned, double*)> f, int z) {
    f(a, domainVariable[x], 1, (N+1)*(N+1)*ne_x*ne_y, domainVariable[z]);

    return;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable z to f(x, y).
 *
 * @Param a Scaling value for x
 * @Param x The first parameter of the function.
 * @Param b Scaling value for y
 * @Param y The second parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForVariables(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, unsigned, double*)> f, int z) {
    f(a, domainVariable[x], b, domainVariable[y], 1, (N+1)*(N+1)*ne_x*ne_y, domainVariable[z]);

    return;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable z to f(w, x, y).
 *
 * @Param a Scaling for first parameter.
 * @Param w The first parameter of the function
 * @Param b Scaling for second parameter.
 * @Param x The second parameter of the function.
 * @Param c Scaling for third parameter.
 * @Param y The third parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForVariables(double a, int w, double b, int x, double c, int y, function<void(double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z) {
    f(a, domainVariable[w], b, domainVariable[x], c, domainVariable[y], 1, (N+1)*(N+1)*ne_x*ne_y, domainVariable[z]);
    
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable z to f(a, b, c, d).
 *
 * @Param t Scaling for first parameter.
 * @Param a The first parameter of the function
 * @Param u Scaling for second parameter.
 * @Param b The second parameter of the function.
 * @Param v Scaling for third parameter.
 * @Param c The third parameter of the function.
 * @Param x Scaling for fourth parameter.
 * @Param d The fourth parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForVariables(double t, int a, double u, int b, double v, int c, double x, int d, function<void(double, double*, double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z) {
    f(t, domainVariable[a], u, domainVariable[b], v, domainVariable[c], x, domainVariable[d], 1, (N+1)*(N+1)*ne_x*ne_y, domainVariable[z]);
    
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable z to f(a, b, c, d, e).
 *
 * @Param p Scaling for first parameter.
 * @Param a The first parameter of the function
 * @Param q Scaling for second parameter.
 * @Param b The second parameter of the function.
 * @Param r Scaling for third parameter.
 * @Param c The third parameter of the function.
 * @Param s Scaling for fourth parameter.
 * @Param d The fourth parameter of the function.
 * @Param t Scaling for fifth parameter.
 * @Param e The fifth parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForVariables(double p, int a, double q, int b, double r, int c, double s, int d, double t, int e, function<void(double, double*, double, double*, double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z){
    f(p, domainVariable[a], q, domainVariable[b], r, domainVariable[c], s, domainVariable[d], t, domainVariable[e], 1, (N+1)*(N+1)*ne_x*ne_y, domainVariable[z]);
    
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable z to f(a, b, c, d, e, g).
 *
 * @Param p Scaling for first parameter.
 * @Param a The first parameter of the function
 * @Param q Scaling for second parameter.
 * @Param b The second parameter of the function.
 * @Param r Scaling for third parameter.
 * @Param c The third parameter of the function.
 * @Param s Scaling for fourth parameter.
 * @Param d The fourth parameter of the function.
 * @Param t Scaling for fifth parameter.
 * @Param e The fifth parameter of the function.
 * @Param u Scaling for sixth parameter.
 * @Param g The sixth paramter.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForVariables(double p, int a, double q, int b, double r, int c, double s, int d, double t, int e, double u, int g, function<void(double, double*, double, double*, double, double*, double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z){
    f(p, domainVariable[a], q, domainVariable[b], r, domainVariable[c], s, domainVariable[d], t, domainVariable[e], u, domainVariable[g], 1, (N+1)*(N+1)*ne_x*ne_y, domainVariable[z]);
    
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable z to f(a*x).
 *
 * @Param a Scaling for first parameter.
 * @Param x The first parameter of the function
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsforDomainVariablesfromCellCenterVariables(double a, int x, function<void(double, double*, unsigned, unsigned, double*, unsigned)> f, int z){
    f(a, cellcenterVariable[x], ne_x*ne_y, (N+1)*(N+1)*ne_x*ne_y, domainVariable[z], (N+1)*(N+1));
    return ;
}
/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of Boundary variable z to f(x, y).
 *
 * @Param a Scaling value for x
 * @Param x The first parameter of the function.
 * @Param b Scaling value for y
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForBoundaryVariables(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, unsigned, double*)> f, int z) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            elements[i][j]->setFunctionsForBoundaryVariables(a, x, b, y, f, z);
    return;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of Boundary variable z to f(w, x, y).
 *
 * @Param a Scaling value for w
 * @Param w The first parameter of the function.
 * @Param b Scaling value for x
 * @Param x The second parameter of the function.
 * @Param c Scaling value for y
 * @Param y The third parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForBoundaryVariables(double a, int w, double b, int x, double c, int y, function<void(double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            elements[i][j]->setFunctionsForBoundaryVariables(a, w, b, x, c, y, f, z);
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of variable stored at center z to f(x, y).
 *
 * @Param x The first parameter of the function.
 * @Param y The second parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 *                  with arguments - a, x, b, y, Sizeof(x), z  
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForVariablesCellCentered(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, double*)> f, int z) {
    f(a, cellcenterVariable[x], b, cellcenterVariable[y], ne_x*ne_y, cellcenterVariable[z]);

    return;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of cell center variable stored at z to f(x, y) dependent on domain variables x and y.
 *
 * @Param x The first parameter of the function.
 * @Param y The second parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 *                  with arguments - a, x, b, y, Sizeof(x), Sizeof(z), z, Nodes in element   
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::setFunctionsForCellCenterVariablesfromDomainVariables(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, unsigned, double*, unsigned)> f, int z) {
    f(a, domainVariable[x], b, domainVariable[y], (N+1)*(N+1)*ne_x*ne_y, ne_x*ne_y, cellcenterVariable[z], (N+1)*(N+1));

    return;
}





double DG_Field_2d::l2Norm(int v1, int v2) {
    double norm = 0.0;
    double elementNorm;
    int j;
    #pragma omp parallel for private(j, elementNorm) reduction(+:norm)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++){
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
 * @Param level This int is used to identify the level of limiting rquired.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::checkPositivity(int v, int cm, string level) {
     int a, k;
     //if( level == "One")
      {
         double *domVar = domainVariable[v]; 
         int j;
         #pragma omp parallel for private(j,a,k)
         for(int i=0; i <ne_x; ++i) {
             for(j=0; j<ne_y; ++j) {
                 a = (N+1)*(N+1)*(i*ne_y + j);
                 k = 0;
                 while( (k < (N+1)*(N+1)) && !PositivityMarker[i*ne_y + j] ) {
                     if(domVar[a + k] < 0) PositivityMarker[i*ne_y + j] = true;
                     k++;
                 }
             }
         }
     }
     /*for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
           // if ( cellcenterVariable[cm][i*ne_y + j]) 
            { 
                if ( level == "One" )
                 {
                    PositivityMarker[i*ne_y + j] = elements[i][j]->checkPositivity(v, level);
                }
                else if(!PositivityMarker[i*ne_y + j] ) {
                    PositivityMarker[i*ne_y + j] = elements[i][j]->checkPositivity(v, level);
                }
        } */
            
    return;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function resets the positivity markers in cell.
 * 
 * @Param v Either true or false
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::resetPositivity( bool v) {
    int j;
    #pragma omp parallel for private(j)
     for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            PositivityMarker[i*ne_y + j] = v;

    return;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function adds Conservative Variables.
 *
 * @Param v The corresponding conservative variable in the system being solved
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::addConservativeVariables(int v) {
    for(int i = 0; i < ne_x; i++)
    for (int j = 0; j < ne_y; j++) {
      elements[i][j]->addConservativeVariables(v);
    }

    return ;
}

void DG_Field_2d::addConservativeVariables(vector<int> V) {
    for(int i = 0; i < ne_x; i++)
    for (int j = 0; j < ne_y; j++) {
      elements[i][j]->addConservativeVariables(V);
    }

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function update Boundary values after every timestep.
 *
 * @Param time the solution time.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::updateBoundary(double time) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
    for (j = 0; j < ne_y; j++) {
      elements[i][j]->updateBoundary(time);
    }

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the maximum value in the Field .
 *
 * @Param v This is a int which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
double DG_Field_2d::FindMax(int v) {
    double *domVar = domainVariable[v];
    double Max = domVar[0]; 
    #pragma omp parallel for reduction(max:Max)
    for(int i=0; i< (N+1)*(N+1)*ne_x*ne_y; ++i) {
        Max  = max(Max, domVar[i]);
    }
    return Max;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the maximum value in the Field for Cell Centered variables .
 *
 * @Param v This is a int which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
double DG_Field_2d::FindMaxCellCentered(int v) {
    double *CCVar = cellcenterVariable[v];
    double Max = CCVar[0]; 
    #pragma omp parallel for reduction(max:Max)
    for(int i=0; i< ne_x*ne_y; ++i) {
        Max  = max(Max, CCVar[i]);
    }
    return Max;
}
/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the minimum dx in the field.
 *
 */
/* ----------------------------------------------------------------------------*/
double DG_Field_2d::FindMindx() {
    double Min = min(elements[0][0]->dxMin, elements[0][0]->dyMin) ; 
    int j;
    #pragma omp parallel for private(j) reduction(min:Min)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            Min = min(min(elements[i][j]->dxMin, elements[i][j]->dyMin), Min);
    return Min;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the minimum dx in the field.
 *
 * @Param dx To store the minimum dx in the element
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::FindMindx(int dx) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            cellcenterVariable[dx][i*ne_y + j] = min(elements[i][j]->dxMin, elements[i][j]->dyMin);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the minimum dt in the field.
 *
 * @Param dt To find th mimimum dt from all dt.
 */
/* ----------------------------------------------------------------------------*/
double DG_Field_2d::FindMindt(int dt) {
    double *CCVar = cellcenterVariable[dt];
    double Mini = CCVar[0];
    int j;
    #pragma omp parallel for private(j) reduction(min:Mini)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            Mini = min(CCVar[i*ne_y + j] , Mini);

    return Mini;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the maximum U in each element int the field.
 *
 * @Param dx To store the max U in the element
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::FindUMax(int U, int V, int D, int T, int UMax) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            cellcenterVariable[UMax][i*ne_y + j] = elements[i][j]->FindUMax(U, V, D, T);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the minimum dt in the field.
 *
 * @Param dx To store the minimum dx in the element
 * @Param dt Stores the minimum dt in an element
 * @Param U Stores the Maximum eigen value in the element
 * @Param CFL CFL to be used
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::FindTimestep(int dt, int dx, int U, double CFL) {
    double *CCVar_dt = cellcenterVariable[dt], *CCVar_dx = cellcenterVariable[dx], *CCVarU = cellcenterVariable[U];
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            //CCVar_dt[i*ne_y + j] = CFL * (CCVar_dx[i*ne_y + j]/ CCVarU[i*ne_y + j] );
            CCVar_dt[i*ne_y + j] = CFL * 1.0/ CCVarU[i*ne_y + j];

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the total runtime of the simulation.
 *
 * @Param t T Final.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::SetFinalTime(double t) {
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            elements[i][j]->SetFinalTime(t);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates the runtime of the simulation.
 *
 * @Param dt Stores the minimum dt in an element
 * @Param Time Current dt
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::updateTime(int dt, int Time) {
    double *CCVar_dt = cellcenterVariable[dt];
    int j;
    #pragma omp parallel for private(j)
    for(int i = 0; i < ne_x; i++)
        for(j = 0; j < ne_y; j++)
            elements[i][j]->updateTime(CCVar_dt[i*ne_y + j], Time);

    return ;
}