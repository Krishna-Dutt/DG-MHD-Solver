#include "../../includes/DG_Element_2d/DG_Element_2d.h"
#include "../../includes/Utilities/Zeros.h"

#include "../../includes/Utilities/Inverse.h"

#include "../../includes/Utilities/FluxMatrix.h"
#include "../../includes/Utilities/MassMatrix.h"
#include "../../includes/Utilities/DerivativeMatrix.h"
#include "../../includes/Utilities/LobattoNodes.h"
#include "../../includes/Utilities/MinMod.h"
#include "../../includes/Utilities/MathOperators.h"
#include "../../includes/Utilities/MaterialProperties.h"

#include <cmath>

#define MAX(a, b)(a>b?a:b)
#define MIN(a, b)(a<b?a:b)

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the constructor which initializes the member variables of the class. It also populates the arrays
 * X, Y as once the co-ordinates of the rectangle is known these arrays can be computed.
 *
 * @Param[in] _N The order of the polynomial which is to be used.
 * @Param[in] x1 The x-coord of the bottom left corner.
 * @Param[in] y1 The y-coord of the bottom left corner.
 * @Param[in] x2 The x-coord of the top right corner.
 * @Param[in] y2 The y-coord of the top right corner.
 */
/* ----------------------------------------------------------------------------*/
DG_Element_2d::DG_Element_2d(int _N, double x1, double y1, double x2, double y2) {
    N = _N;
    x_start = x1;
    y_start = y1;
    x_end   = x2;
    y_end   = y2;


    X = new double[(N+1)*(N+1)];
    Y = new double[(N+1)*(N+1)];

    double* nodes = new double[N+1]; // Alloting space for the lobatto nodes.
    
    lobattoNodes(nodes, N+1); // Found the lobatto nodes in the range -1 to +1.

    for(int k1=0;k1<=N;k1++)
    {
        for(int k2=0;k2<=N;k2++)
        {
            X[k1*(N+1)+k2] = 0.5*(x1 + x2) + 0.5*(x2 - x1)*nodes[k2];
            Y[k1*(N+1)+k2] = 0.5*(y1 + y2) + 0.5*(y2 - y1)*nodes[k1];
        }
    }

    dxMin = abs(X[0]-X[1]);
    dyMin = abs(Y[0]-Y[N+1]);

    delete[] nodes;

    massMatrix = NULL;

    derivativeMatrix_x  =   NULL;
    derivativeMatrix_y  =   NULL; 
    fluxMatrix_top      =   NULL;
    fluxMatrix_right    =   NULL;
    fluxMatrix_bottom   =   NULL;
    fluxMatrix_left     =   NULL;
    inverseMassMatrix   =   NULL;
    topNeighbor         =   this;
    rightNeighbor       =   this;
    leftNeighbor        =   this;
    bottomNeighbor      =   this;

    vanderMandMatrix    =   NULL;
    inverseVanderMandMatrix = NULL;

    PositivityMarker   = true;

    System = "EULER";

    //OutFlow.resize(0);
    // For Quad element : 0 - Bottom, 1 - Right, 2 - Top, 3 - Left Boundary
    /*OutFlow.push_back(true);
    OutFlow.push_back(true);
    OutFlow.push_back(true);
    OutFlow.push_back(true);
    */
    OutFlow[0] = OutFlow[1] = OutFlow[2] = OutFlow[3] = true;

    RightEigenMatrix = NULL;
    LeftEigenMatrix  = NULL;
    
    Dimension = 4;


}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the destructor which deallocates the member variables of the class and destroys the object.
 *
*/
/* ----------------------------------------------------------------------------*/
DG_Element_2d::~DG_Element_2d() {
    delete[] X;
    delete[] Y;

    /*for ( vector<double*>::iterator itr = variable.begin() ;itr != variable.end(); itr++){
      delete[] (itr->second);
    }
    // Do I need to delete other maps pointing to Boundary elements, since no new memory is dynamically allocated for them ??

    for ( vector<double**>::iterator itr = boundaryTop.begin() ;itr != boundaryTop.end(); itr++){
      delete[] (itr->second);
    }
    for ( vector<double**>::iterator itr = boundaryBottom.begin() ;itr != boundaryBottom.end(); itr++){
      delete[] (itr->second);
    }
    for ( vector<double**>::iterator itr = boundaryRight.begin() ;itr != boundaryRight.end(); itr++){
      delete[] (itr->second);
    }
    for ( vector<double**>::iterator itr = boundaryLeft.begin() ;itr != boundaryLeft.end(); itr++){
      delete[] (itr->second);
    }

    for ( vector<double**>::iterator itr = neighboringTop.begin() ;itr != neighboringTop.end(); itr++){
      delete[] (itr->second);
    }
    for ( vector<double**>::iterator itr = neighboringBottom.begin() ;itr != neighboringBottom.end(); itr++){
      delete[] (itr->second);
    }
    for ( vector<double**>::iterator itr = neighboringRight.begin() ;itr != neighboringRight.end(); itr++){
      delete[] (itr->second);
    }
    for ( vector<double**>::iterator itr = neighboringLeft.begin() ;itr != neighboringLeft.end(); itr++){
      delete[] (itr->second);
    }*/
}



/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function which deallocates the matrices of the class.
 *
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::Destroy_Matrices() {
    delete[] massMatrix;
    delete[] derivativeMatrix_x;
    delete[] derivativeMatrix_y;
    delete[] transposeDerivateMatrix_x;
    delete[] transposeDerivateMatrix_y;
    delete[] fluxMatrix_top;
    delete[] fluxMatrix_right;
    delete[] fluxMatrix_bottom;
    delete[] fluxMatrix_left;
    delete[] inverseMassMatrix;
    if ( vanderMandMatrix != NULL){
      delete[] vanderMandMatrix;
      delete[] inverseVanderMandMatrix;
    }
    if (RightEigenMatrix != NULL){
        delete[] RightEigenMatrix;
        delete[] LeftEigenMatrix;
    }

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This functions adds the name of the governing equation being solved
 *
 * @Param S  This is the name of system being solved.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setSystem(string S) {
    System = S;

    return ;
}
  

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This functions creates space in order to take in one more variable on which operators are needed to be
 * applied.
 *
 * @Param v  This is the name of the variable which is to be added.
 * @Param p  This is the pointer to the allocated memory
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::addVariable_withoutBoundary(int v, double *p) {
    //double * newVariable = new double[(N+1)*(N+1)]; /// Allocating the space for the new variable which is to be created.
    /// Allocating the space for the new variable which is to be created.
    if (v == 0) variable.clear();

    variable.push_back(p); /// Now assigning the same to the map.

    double *b_top       = NULL; 
    double *b_bottom    = NULL;
    double *b_left      = NULL;
    double *b_right     = NULL;

   
    double *n_top       = NULL;
    double *n_left      = NULL;
    double *n_right     = NULL;
    double *n_bottom    = NULL;

    if (v == 0) {
        boundaryTop.clear();
        boundaryBottom.clear();
        boundaryLeft.clear();
        boundaryRight.clear();

        neighboringTop.clear();
        neighboringBottom.clear();
        neighboringLeft.clear();
        neighboringRight.clear();
    } 

    boundaryTop.push_back(b_top);
    boundaryRight.push_back(b_right);
    boundaryBottom.push_back(b_bottom);
    boundaryLeft.push_back(b_left);


    neighboringTop.push_back(n_top);
    neighboringRight.push_back(n_right);
    neighboringBottom.push_back(n_bottom);
    neighboringLeft.push_back(n_left);


    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This functions creates space in order to take in one more variable on which operators are needed to be
 * applied.
 *
 * @Param v  This is the name of the variable which is to be added.
 * @Param p  Pointer to memory allocation for variable.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::addVariable_withBoundary(int v, double *p) {
    if (v == 0) variable.clear();

    variable.push_back(p); /// Now assigning the same to the map.
    
    // **b_top is used because it will store the address to the boundary element, So whenever the actual value in the double* of the variable is changed then this will also change automatically. The same holds for all the other following mentioned variables.
    
    double *b_top    ; 
    double *b_bottom ;
    double *b_left   ;
    double *b_right  ;

   
    double *n_top    ;
    double *n_left   ;
    double *n_right  ;
    double *n_bottom ;

    
        n_bottom =   b_bottom = &(variable[v][0]);    
        n_top    =   b_top    = &(variable[v][N*(N+1)+0]);
        n_left   =   b_left   = &(variable[v][0*(N+1)]);
        n_right  =   b_right  = &(variable[v][0*(N+1)+N]);
    

    if (v == 0) {
        boundaryTop.clear();
        boundaryBottom.clear();
        boundaryLeft.clear();
        boundaryRight.clear();

        neighboringTop.clear();
        neighboringBottom.clear();
        neighboringLeft.clear();
        neighboringRight.clear();
    } 

    boundaryTop.push_back(b_top);
    boundaryRight.push_back(b_right);
    boundaryBottom.push_back(b_bottom);
    boundaryLeft.push_back(b_left);


    neighboringTop.push_back(n_top);
    neighboringRight.push_back(n_right);
    neighboringBottom.push_back(n_bottom);
    neighboringLeft.push_back(n_left);

    boundaryVariables.push_back(v);

    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This functions creates space in order to a cell centered variable on which operators are needed to be
 * applied.
 *
 * @Param v  This is the name of the variable which is to be added.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::addVariable_CellCentered(int v) {
    double * newVariable = new double[1] ;/// Allocating the space for the new variable which is to be created.
    variable[v] = newVariable; /// Now assigning the same to the map.

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This functions resets cell centered variable on which operators are needed to be
 * applied.
 *
 * @Param v  This is the name of the variable which is to be reset.
 * @Param value The value to which the variable is to be reset.
 */
/* ----------------------------------------------------------------------------*/
/*void DG_Element_2d::ResetVariables_CellCentered(int v, double value) {
   *variable[v] = value; 
    return ;
}*/

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function resets the Map to the Outflow Boundaries.
 *
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::ResetMap_OutFlow() {
    for(int i=0; i<4; ++i) {
        OutFlow[i] = false;
    }
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates the Map to the Outflow Boundaries.
 *
 * @Param u This is the velocity in the x direction.
 * @Param v This is the velocity in the y direction.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::updateOutFlowBoundary(int u, int v) {
  double start = -1.0;
  double end = 1.0;
  
  if ( Sum(1, N+1, boundaryTop[v]) < 0.0) {
    OutFlow[2] = true;
  }
  if ( Sum(1, N+1, boundaryBottom[v]) > 0.0) {
    OutFlow[0] = true;
  }
  
  if ( Sum(N+1, (N+1)*(N+1), boundaryLeft[v]) > 0.0) {
    OutFlow[3] = true;
  }
  if (Sum(N+1, (N+1)*(N+1)-N, boundaryRight[v]) < 0.0) {
    OutFlow[1] = true;
  }

  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates cell markers for Shock Detection, using LXRCF.
 *
 * @Param v This is the quantity used for detection.
 * @Param m This is the variables used to store the cell marker.
*/
/* ----------------------------------------------------------------------------*/
double DG_Element_2d::updateCellMarker(int v) {
  double radius = 1.0;
  double OutflowSize = 0.0;
  double MaxVariable = 0.0;
  double VariableFlux = 0.0;
  double Marker = 0.0;

  double *NumFlux = new double[N+1];
  
  zeros(NumFlux, N+1);
  if (OutFlow[2]) {
      Addabc(1.0, boundaryTop[v], -1.0, neighboringTop[v], 1.0, NumFlux, 1, N+1, NumFlux);
      OutflowSize += abs(x_end -x_start);
  }
  if (OutFlow[0]) {
      Addabc(1.0, boundaryBottom[v], -1.0, neighboringBottom[v], 1.0, NumFlux, 1, N+1, NumFlux);
      OutflowSize += abs(x_end -x_start);
  }
  if (OutFlow[0] || OutFlow[2]) {
      VariableFlux += lobattoIntegration(x_start, x_end, N, 1, NumFlux);
  }

  zeros(NumFlux, N+1);
  if (OutFlow[3]) {
      for(int i=0; i<N+1; ++i) {
          NumFlux[i] += boundaryLeft[v][i*(N+1)]; 
      }
      for(int i=0; i<N+1; ++i) {
          NumFlux[i] -= neighboringLeft[v][i*(N+1)]; 
      }
      OutflowSize += abs(y_end -y_start);
  }
  if (OutFlow[1]) {
      for(int i=0; i<N+1; ++i) {
          NumFlux[i] += boundaryRight[v][i*(N+1)]; 
      }
      for(int i=0; i<N+1; ++i) {
          NumFlux[i] -= neighboringRight[v][i*(N+1)]; 
      }
      OutflowSize += abs(y_end -y_start);
  }
  if(OutFlow[1] || OutFlow[3]) {
      VariableFlux += lobattoIntegration(y_start, y_end, N, 1, NumFlux);
  }

  /*if (OutFlow[2]) {
    VariableFlux += lobattoIntegration(x_start, x_end, N, 1, boundaryTop[v]);
    VariableFlux -= lobattoIntegration(x_start, x_end, N, 1, neighboringTop[v]);
    OutflowSize += abs(x_end -x_start);

  }
  if (OutFlow[0]) {
    VariableFlux += lobattoIntegration(x_start, x_end, N, 1, boundaryBottom[v]);
    VariableFlux -= lobattoIntegration(x_start, x_end, N, 1, neighboringBottom[v]);
    OutflowSize += abs(x_end -x_start);

  }
  if (OutFlow[3]) {
    VariableFlux += lobattoIntegration(y_start, y_end, N, N+1, boundaryLeft[v]);
    VariableFlux -= lobattoIntegration(y_start, y_end, N, N+1, neighboringLeft[v]);
    OutflowSize += abs(y_end -y_start);

  }
  if (OutFlow[1]) {
    VariableFlux += lobattoIntegration(y_start, y_end, N, N+1, boundaryRight[v]);
    VariableFlux -= lobattoIntegration(y_start, y_end, N, N+1, neighboringRight[v]);
    OutflowSize += abs(y_end -y_start);

  }*/
  

  MaxVariable = variable[v][0];
  for (int i = 0; i < (N+1)*(N+1) ; ++i) {
    MaxVariable = MAX(MaxVariable,variable[v][i]);
  }


// Assert that MaxVariable is never equal to ZERO !!  

  radius = MIN(abs(x_start-x_end),abs(y_start-y_end)) * 0.5;

  Marker = abs(VariableFlux) / ( abs(OutflowSize) * MaxVariable * pow(radius, 0.5 * (N+1)));
  if (Marker > 1.0) {
    Marker = 1.0;
  }
  else {
    Marker = 0.0;
  }
  
  for(int b=0; b <(N+1)*(N+1); ++b) {
      variable[35][b] = Marker; // CellMarkerG
  }
 
 delete[] NumFlux;

 return Marker;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This functions checks the positivity of data stored in cell.
 *
 * @Param v  This is the name of the variable whose positivity is to be checked.
 * @Param cm This is the cell marker used to identify troubled cells.
 *v@Param level This string is used to identify the level of limiting required.
 */
/* ----------------------------------------------------------------------------*/
bool DG_Element_2d::checkPositivity(int v, string level) {     
        for (int i=0; i< (N+1)*(N+1); ++i) {
            if (variable[v][i] < 0.0)
             {
                return true ;
            }
        }

       /* if (level == "One"){
            //PositivityMarker = false;
            return false;
        }
        */         
    return false;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function resets the positivity markers in cell.
 *
 */
/* ----------------------------------------------------------------------------*/
/*void DG_Element_2d::resetPositivity() {
    PositivityMarker = true;
         
    return;
}*/


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function computes the moments for a given variable using VandeMandMatrix with Legendre Basis.
 *
 * @Param v This is the variable whose moments are to be captured.
 * @Param m This is the variable used to store the computed moments.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::computeMoments(int v, int m) {
  /// Multiplying inverse of VanderMand Matrix with the variable array to get the corresponding moments.
  cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 1.0, inverseVanderMandMatrix,(N+1)*(N+1), variable[v],1,0,variable[m],1);
  
  return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function computes the nodals values of the variable given its moments
 * using VandeMandMatrix with Legendre Basis.
 *
 * @Param m This gives the moments of the variable
 * @Param v This is the variable whose nodal values are to be computed.
 * @Param cm This is the cell marker used to identify troubled cells.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::convertMomentToVariable(int m, int v) {
  /// Multiplying  VanderMand Matrix with the moments to obtained the nodal values of the variable.

  cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 1.0, vanderMandMatrix,(N+1)*(N+1), variable[m],1,0,variable[v],1);

  return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function limits the moments 
 * using Lilia's Moment Limiter.
 *
 * @Param m This gives the moments of the variable
 * @Param modm This is the variable to store modified moments.
 * @Param cm This is the cell marker used to identified troubled cells.
 * @Param Index This is the index to start the limiting process.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::limitMoments(int m, int modm, unsigned Index) {

  { // Checking if cell marker is not equal to zero
    int count, Tempi, Tempj, i, j;
    count = N+1;
    double Temp1, Temp2, AlphaN;
    double epsilon = 1e-13;
    AlphaN = sqrt((2.0*N -1.0)/(2.0*N +1));  // Similar to a diffusion coefficient
    /*vector<double> Var1, Var2 ;
    double M = 50.0;
    double min_dx = min(abs(X[0]-X[1]), abs(Y[0]-Y[1]));
    Var1.push_back(M*min_dx*min_dx);
    Var2.push_back(M*min_dx*min_dx);
    */  
    
    // Ensuring that Cell avergae remains the  same after limiting !!
    variable[modm][0] = variable[m][0];

    if ( N == 1) {
        epsilon = 0.0 ;
    }
    else {
        epsilon = 1e-16;
    }

    for(i=Index; i > 0; i = i - (N+2)) {
     --count;
     //AlphaN = sqrt((2.0*(count)-1.0)/(2.0*(count)+1.0));
     for(j=0; j < count; ++j) {
       Tempi = i-j;
       Tempj = i - j*(N+1);

       /*Var1.resize(0);
       Var2.resize(0);
       Var1.push_back(variable[m][Tempi]);
       Var2.push_back(variable[m][Tempj]);

       //if (abs(rightNeighbor->variable[m][Tempi]) > epsilon ) 
       Var1.push_back(AlphaN*(rightNeighbor->variable[m][Tempi-1] -variable[m][Tempi-1]));
       //if (abs(leftNeighbor->variable[m][Tempi]) > epsilon) 
       Var1.push_back(AlphaN*(variable[m][Tempi-1] -leftNeighbor->variable[m][Tempi-1]));
       //if (abs(topNeighbor->variable[m][Tempi]) > epsilon ) 
       Var1.push_back(AlphaN*(topNeighbor->variable[m][Tempi-(N+1)] -variable[m][Tempi-(N+1)]));
       //if (abs(bottomNeighbor->variable[m][Tempi]) > epsilon) 
       Var1.push_back(AlphaN*(variable[m][Tempi-(N+1)] -bottomNeighbor->variable[m][Tempi-(N+1)])); 
       
       //if (abs(rightNeighbor->variable[m][Tempj]) > epsilon ) 
       Var2.push_back(AlphaN*(rightNeighbor->variable[m][Tempj-1] -variable[m][Tempj-1]));
       //if (abs(leftNeighbor->variable[m][Tempj]) > epsilon) 
       Var2.push_back(AlphaN*(variable[m][Tempj-1] -leftNeighbor->variable[m][Tempj-1]));
       //if (abs(topNeighbor->variable[m][Tempj]) > epsilon ) 
       Var2.push_back(AlphaN*(topNeighbor->variable[m][Tempj-(N+1)] -variable[m][Tempj-(N+1)]));
       //if (abs(bottomNeighbor->variable[m][Tempj]) > epsilon) 
       Var2.push_back(AlphaN*(variable[m][Tempj-(N+1)] -bottomNeighbor->variable[m][Tempj-(N+1)])); 
       
       Temp1 = MinMod(Var1);
       Temp2 = MinMod(Var2);
       */

       // Original minmod detector
       Temp1 = MinMod(variable[m][Tempi], AlphaN*(rightNeighbor->variable[m][Tempi-1] -variable[m][Tempi-1]), AlphaN*(variable[m][Tempi-1] -leftNeighbor->variable[m][Tempi-1]) , AlphaN*(topNeighbor->variable[m][Tempi-(N+1)] -variable[m][Tempi-(N+1)]), AlphaN*(variable[m][Tempi-(N+1)] -bottomNeighbor->variable[m][Tempi-(N+1)]));
       Temp2 = MinMod(variable[m][Tempj], AlphaN*(rightNeighbor->variable[m][Tempj-1] -variable[m][Tempj-1]), AlphaN*(variable[m][Tempj-1] -leftNeighbor->variable[m][Tempj-1]) , AlphaN*(topNeighbor->variable[m][Tempj-(N+1)] -variable[m][Tempj-(N+1)]), AlphaN*(variable[m][Tempj-(N+1)] -bottomNeighbor->variable[m][Tempj-(N+1)]));
       
       if (abs(Temp1-variable[modm][Tempi]) > epsilon || abs(Temp2-variable[modm][Tempj]) > epsilon ) {
         variable[modm][Tempi] = Temp1;
         variable[modm][Tempj] = Temp2;
       }
       else 
       {
         return ; // Need to exit both loops
       }
     }
     // Special Case for end values, when Tempi or Tempj access zero order polynomials !!
       Tempi = i-j;
       Tempj = i - j*(N+1);


       /*Var1.resize(0);
       Var2.resize(0);
       Var1.push_back(variable[m][Tempi]);
       Var2.push_back(variable[m][Tempj]);

       //if (abs(topNeighbor->variable[m][Tempi]) > epsilon ) 
       Var1.push_back(AlphaN*(topNeighbor->variable[m][Tempi-(N+1)] -variable[m][Tempi-(N+1)]));
       //if (abs(bottomNeighbor->variable[m][Tempi]) > epsilon) 
       Var1.push_back(AlphaN*(variable[m][Tempi-(N+1)] -bottomNeighbor->variable[m][Tempi-(N+1)])); 
       
       //if (abs(rightNeighbor->variable[m][Tempj]) > epsilon ) 
       Var2.push_back(AlphaN*(rightNeighbor->variable[m][Tempj-1] -variable[m][Tempj-1]));
       //if (abs(leftNeighbor->variable[m][Tempj]) > epsilon) 
       Var2.push_back(AlphaN*(variable[m][Tempj-1] -leftNeighbor->variable[m][Tempj-1]));
       
       Temp1 = MinMod(Var1);
       Temp2 = MinMod(Var2);
       */
       // Original Detector
       Temp1 = MinMod(variable[m][Tempi], AlphaN*(topNeighbor->variable[m][Tempi-(N+1)] -variable[m][Tempi-(N+1)]), AlphaN*(variable[m][Tempi-(N+1)] -bottomNeighbor->variable[m][Tempi-(N+1)]));
       Temp2 = MinMod(variable[m][Tempj], AlphaN*(rightNeighbor->variable[m][Tempj-1] -variable[m][Tempj-1]), AlphaN*(variable[m][Tempj-1] -leftNeighbor->variable[m][Tempj-1]));
       
       if ( abs(Temp1-variable[modm][Tempi]) > epsilon || abs(Temp2-variable[modm][Tempj]) > epsilon ) {
         variable[modm][Tempi] = Temp1;
         variable[modm][Tempj] = Temp2;
       }
       else //if( Temp1 !=0 && Temp2 !=0)
       {
         return ; // Need to exit both loops
       }
    }

 }

  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to limit the Moments for given variable, using Characteristic Limiter.
 * 
 * @Param array V array of moments of Conservative Variables.
 * @Param C   moment of characteristic variables, to store modified moments.
 * @Param Index Index correspoding to Characteristic variable.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::limitMoments(int *V, int C, unsigned Index) {
  { 
    int count, Tempi, Tempj, i, j;
    double Temp1, Temp2, AlphaN;
    double epsilon = 1e-13;
    count = N+1;
    AlphaN = sqrt((2.0*N -1.0)/(2.0*N +1));  // Similar to a diffusion coefficient
    vector<double> Var1, Var2 ;
    double Sum1, Sum2, Sum3, Sum4;
    
    if ( N == 1) {
        epsilon = 0.0 ;
    }
    else {
        epsilon = 1e-16;
    }

        for(i=(N+1)*(N+1)-1; i > 0; i = i - (N+2)) {
          --count;
          AlphaN = sqrt((2.0*(count)-1.0)/(2.0*(count)+1.0));
          //AlphaN = 0.5/sqrt(4.0*count*count -1.0);
          //AlphaN = 0.25*(4.0*count-1.0)/sqrt(4.0*count*count -1.0);
          for(j=0; j < count; ++j) {
             Tempi = i-j;
             Tempj = i - j*(N+1);

             Sum1 = Sum2 = Sum3 = Sum4 = 0;
             for(int z=0; z<Dimension; ++z) {
                 Sum1 += LeftEigenMatrix[Index*Dimension+z] * rightNeighbor->variable[V[z]][Tempi-1];
                 Sum2 += LeftEigenMatrix[Index*Dimension+z] * leftNeighbor->variable[V[z]][Tempi-1]; 
                 Sum3 += LeftEigenMatrix[Index*Dimension+z] * topNeighbor->variable[V[z]][Tempi-(N+1)];
                 Sum4 += LeftEigenMatrix[Index*Dimension+z] * bottomNeighbor->variable[V[z]][Tempi-(N+1)];   
             }

             Temp1 = MinMod(variable[C][Tempi], AlphaN*(Sum1 -variable[C][Tempi-1]), AlphaN*(variable[C][Tempi-1] -Sum2), AlphaN*(Sum3 -variable[C][Tempi-(N+1)]), AlphaN*(variable[C][Tempi-(N+1)] -Sum4));
                          
             Sum1 = Sum2 = Sum3 = Sum4 = 0;
             for(int z=0; z<Dimension; ++z) {
                 Sum1 += LeftEigenMatrix[Index*Dimension+z] * rightNeighbor->variable[V[z]][Tempj-1];
                 Sum2 += LeftEigenMatrix[Index*Dimension+z] * leftNeighbor->variable[V[z]][Tempj-1]; 
                 Sum3 += LeftEigenMatrix[Index*Dimension+z] * topNeighbor->variable[V[z]][Tempj-(N+1)];
                 Sum4 += LeftEigenMatrix[Index*Dimension+z] * bottomNeighbor->variable[V[z]][Tempj-(N+1)]; 
             }

             Temp2 = MinMod(variable[C][Tempj], AlphaN*(Sum1 -variable[C][Tempj-1]), AlphaN*(variable[C][Tempj-1] -Sum2), AlphaN*(Sum3 -variable[C][Tempj-(N+1)]), AlphaN*(variable[C][Tempj-(N+1)] -Sum4));

             if (abs(Temp1-variable[C][Tempi]) > epsilon || abs(Temp2-variable[C][Tempj]) > epsilon ) {
                 variable[C][Tempi] = Temp1;
                 variable[C][Tempj] = Temp2;
             }
             else {
                 return ;
             }
          } 
             // Special Case for end values, when Tempi or Tempj access zero order polynomials !!
             Tempi = i-j;
             Tempj = i - j*(N+1);
    
             Sum1 = Sum2 = Sum3 = Sum4 = 0;
             for(int z=0; z<Dimension; ++z) {
                 Sum1 += LeftEigenMatrix[Index*Dimension+z] * topNeighbor->variable[V[z]][Tempi-(N+1)]; 
                 Sum2 += LeftEigenMatrix[Index*Dimension+z] * bottomNeighbor->variable[V[z]][Tempi-(N+1)]; 
                 Sum3 += LeftEigenMatrix[Index*Dimension+z] * rightNeighbor->variable[V[z]][Tempj-1]; 
                 Sum4 += LeftEigenMatrix[Index*Dimension+z] * leftNeighbor->variable[V[z]][Tempj-1]; 
             }
            
             Temp1 = MinMod(variable[C][Tempi], AlphaN*(Sum1 -variable[C][Tempi-(N+1)]), AlphaN*(variable[C][Tempi-(N+1)] -Sum2));
             Temp2 = MinMod(variable[C][Tempj], AlphaN*(Sum3 -variable[C][Tempj-1]), AlphaN*(variable[C][Tempj-1] -Sum4));
             
             if ( abs(Temp1-variable[C][Tempi]) > epsilon || abs(Temp2-variable[C][Tempj]) > epsilon ) {
                variable[C][Tempi] = Temp1;
                variable[C][Tempj] = Temp2;
             }
             else {
                 return ;
             }
        }

   }

  return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function which sets the Right and Left Eigen Matrix corresponding to directionalong flow
 * 
 * @Param dimension No of systme of Equations
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setEigenMatrices(unsigned dimension) {
  
 RightEigenMatrix = new double[dimension*dimension];
 LeftEigenMatrix = new double[dimension*dimension];

 Dimension = dimension;

  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to find  the Right and Left Eigen Matrix corresponding to direction along flow
 * 
 * @Param array V array of index corresponding to field variables required to compute eigen matrices.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::findEigenMatrices(int *V) {
  double u, v, c, H, P;
  double epsilon = 1e-10;
  double q, qn, nx, ny, qe;
  double DVx, DVy, D, DE, dPdx, dPdy;
  DVx = variable[V[0]][0];
  DVy = variable[V[1]][0];
  D   = variable[V[2]][0];
  DE  = variable[V[3]][0];
  dPdx = variable[V[4]][0];
  dPdy = variable[V[5]][0];
  
  u = DVx/D;
  v = DVy/D;
  //dPdx = variable["dPdxMoment"][0];
  //dPdy = variable["dPdyMoment"][0];

 /* if ( abs(u*v) > epsilon || abs(u) > epsilon || abs(v) > epsilon ) {
      nx = u/sqrt(u*u + v*v);
      ny = v/sqrt(u*u + v*v);
  }
  else 
  {
      nx = 1.0; // Asuming 1D flow in x direction
      ny = 0.0;
  }*/
  if (  abs(dPdx) > epsilon || abs(dPdy) > epsilon ) {
      nx = dPdx/sqrt(dPdx*dPdx + dPdy*dPdy);
      ny = dPdy/sqrt(dPdx*dPdx + dPdy*dPdy);
  }
  else 
  {
      nx = 1.0; // Asuming 1D flow in x direction
      ny = 0.0;
  }

  qn = u*nx + v*ny ;
  q = sqrt(u*u + v*v);
  qe = -u*ny + v*nx;
  P = (gamma-1.0) * ( DE - 0.5*q*q*D ); // Twice pressure
  c = sqrt(gamma*P/D);
  H =  c*c/(gamma-1.0) + 0.5*q*q;  // 0.5 * variable[V[4]][0];// Dividing by 2 to get actual cell average , as normalised Legendre Basis is used to compute moments
  // Change later if this works

  

  // Setting Right Eigen Matrix
  RightEigenMatrix[0] = RightEigenMatrix[1] = RightEigenMatrix[2] = 1.0 ;
  RightEigenMatrix[3] = 0.0;
  RightEigenMatrix[4] = u - c*nx;
  RightEigenMatrix[5] = u;
  RightEigenMatrix[6] = u + c*nx;
  RightEigenMatrix[7] = -ny;
  RightEigenMatrix[8] = v - c*ny;
  RightEigenMatrix[9] = v;
  RightEigenMatrix[10] = v + c*ny;
  RightEigenMatrix[11] = nx;
  RightEigenMatrix[12] = H - qn*c;
  RightEigenMatrix[13] = 0.5 *q*q;
  RightEigenMatrix[14] = H + qn*c;
  RightEigenMatrix[15] = qe;

  /*for(int i=0; i < 4; ++i){
      for(int j=0; j<4; ++j) 
      cout << RightEigenMatrix[i*4+j] << " ";
      cout << "\n";
  }*/

  // Setting Left Eigen Matrix;
  //inverse(RightEigenMatrix,LeftEigenMatrix,Dimension*Dimension);
  LeftEigenMatrix[0] = 0.5 * (0.5*(gamma-1.0)*pow(q/c,2.0) + qn/c);
  LeftEigenMatrix[1] = -0.5*( (gamma -1.0)*u/(c*c) + nx/c);
  LeftEigenMatrix[2] = -0.5*( (gamma -1.0)*v/(c*c) + ny/c);
  LeftEigenMatrix[3] = 0.5 * (gamma-1.0)/(c*c);
  LeftEigenMatrix[4] = 1.0 - 0.5 * (gamma-1.0)*pow(q/c,2.0);
  LeftEigenMatrix[5] = (gamma-1.0)*u/(c*c);
  LeftEigenMatrix[6] = (gamma-1.0)*v/(c*c);
  LeftEigenMatrix[7] = -(gamma-1.0)/(c*c);
  LeftEigenMatrix[8] = 0.5 * (0.5*(gamma-1.0)*pow(q/c,2.0) - qn/c);
  LeftEigenMatrix[9] = -0.5*( (gamma -1.0)*u/(c*c) - nx/c);
  LeftEigenMatrix[10] = -0.5*( (gamma -1.0)*v/(c*c) - ny/c);
  LeftEigenMatrix[11] = 0.5 * (gamma-1.0)/(c*c);
  LeftEigenMatrix[12] = -qe;
  LeftEigenMatrix[13] = -ny;
  LeftEigenMatrix[14] = nx;
  LeftEigenMatrix[15] = 0.0;
  
  /*cout << "\n";
  for(int i=0; i < 4; ++i){
      for(int j=0; j<4; ++j) 
      cout << LeftEigenMatrix[i*4+j] << " ";
      cout << "\n";
  }*/
 
  return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to compute Characteristic Variables.
 * 
 * @Param int V Set of conservative variables..
 * @Param c The Characteristic Variable.
 * @Param I Identifier for Characteristic Variable.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::convertVariabletoCharacteristic(int *V, int c, unsigned I) {
    double Sum = 0.0;
    for(int i=0; i< (N+1)*(N+1); ++i) {
        Sum = 0.0;
        for(int j=0; j < Dimension; ++j){
            Sum += variable[V[j]][i]*LeftEigenMatrix[I+j];
        }     
        variable[c][i] = Sum;
    } 

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to compute Conservative  Variables from Characteristic Variables.
 * 
 * @Param array C Set of characteristic variables..
 * @Param v The Conservative Variable.
 * @Param I Identifier for conservative Variable.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::convertCharacteristictoVariable(int *C, int v, unsigned I) {
    double Sum = 0.0;
    for(int i=0; i< (N+1)*(N+1); ++i) {
        Sum = 0.0;
        for(int j=0; j < Dimension; ++j){
            Sum += variable[C[j]][i]*RightEigenMatrix[I+j];
        }     
        variable[v][i] = Sum;
    }  

    return ;
}



/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This functions creates space in order to take in one more variable, needed only at the Boundary of the cell, on which operators are needed to be
 * applied.
 *
 * @Param v  This is the name of the variable which is to be added.
 */
/* ----------------------------------------------------------------------------*/
/*void DG_Element_2d::addVariable_onlyBoundary(int v) {
    double * newVariable = new double[(N+1)*(4)]; /// Allocating the space for the new variable which is to be created, and stored only at boundary of quadrilateral cell.
    variable[v] = newVariable; /// Now assigning the same to the map.
    
    // **b_top is used because it will store the address to the boundary element, So whenever the actual value in the double* of the variable is changed then this will also change automatically. The same holds for all the other following mentioned variables.
    // Assuming the variable is stored in anticlockwise direction in the 1D array.

    double **b_top       = new double* [N+1];
    double **b_bottom    = new double* [N+1];
    double **b_left      = new double* [N+1];
    double **b_right     = new double* [N+1];

   
    double **n_top       = new double* [N+1];
    double **n_left      = new double* [N+1];
    double **n_right     = new double* [N+1];
    double **n_bottom    = new double* [N+1];

    for(int i=0; i<=N; i++){
        n_bottom[i] =   b_bottom[i] = &(variable[v][i]);
        n_top[i]    =   b_top[i]    = &(variable[v][3*(N+1) -1-i]);
        n_left[i]   =   b_left[i]   = &(variable[v][4*(N+1)-1 -i]);
        n_right[i]  =   b_right[i]  = &(variable[v][i + N+1]);
    }


    boundaryTop[v]      = b_top;
    boundaryRight[v]    = b_right;
    boundaryBottom[v]   = b_bottom;
    boundaryLeft[v]     = b_left;


    neighboringTop[v]   = n_top;
    neighboringRight[v] = n_right;
    neighboringBottom[v]= n_bottom;
    neighboringLeft[v]  = n_left;

    variableOnlyAtBoundary.push_back(v);

    return ;
}*/


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to set the value of a variable with the help of a function.
 *
 * @Param v This is the name of the variable whose value is to be  set.
 * @Param f This is the f(x, y) which is used to initialize the variable.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::initializeVariable(int v, function<double(double, double)> f) {
    for(int i=0; i<((N+1)*(N+1)); i++)
        variable[v][i] = f(X[i], Y[i]);
    
    return ;
}


void DG_Element_2d::setVariableNeighbors(int v) {
    int j;
    if(topNeighbor!=this) {
            neighboringTop[v] = topNeighbor->boundaryBottom[v];
    }
    if(rightNeighbor!=this) {
            neighboringRight[v] = rightNeighbor->boundaryLeft[v];
    }
    if(leftNeighbor!=this) {
            neighboringLeft[v] = leftNeighbor->boundaryRight[v];
    }
    if(bottomNeighbor!=this) {
            neighboringBottom[v] = bottomNeighbor->boundaryTop[v];
    }
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary variables for an element. This must be called whenever a new variable with
 * boundary is added. 
 *
 * @Param v The variable whose information about the boundaries is to be stored
 * @Param type The type of the neighbor whose information is to be looked upon.
 * @Param neighbor The pointer to the element 
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setNeighboringElement(char type, DG_Element_2d* neighbor) {
    switch(type) {
        case 't' : // `t` or `T` for top 
        case 'T' :
            topNeighbor = neighbor;
            break;
        case 'r' : // `r` or `R` for right
        case 'R' :
            rightNeighbor = neighbor;
            break;
        case 'b' : // `b` or `B` for bottom
        case 'B' :
            bottomNeighbor = neighbor;
            break;
        case 'l' : // `l` or `L` for left
        case 'L' :
            leftNeighbor = neighbor;
            break;
        default:
            cout << "WARNING!. No such neighbor type " << type << endl;
    }
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to get the variable `v` differentiated partially w.r.t. `x` and then store it in the
 * variable `vDash`. The function also takes `fluxType` as an input which would describe the numerical scheme that
 * should be used in order to obtain the derivative.
 *
 * @Param v         Variable which is to be differentiated.
 * @Param vDash     Variable in which the derivative is to be stored.
 * @Param conserVar Corresponding Conservative variable.
 * @Param fluxType  The type of flux that is to be used. eg "central"
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::delByDelX(int v, int vDash, int conserVar, string fluxType, int fluxVariable = 999) {
    double dy = (y_end - y_start);
    double dx = (x_end - x_start);
    
    if(fluxType == "central") {
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        zeros(numericalFlux, (N+1)*(N+1));                                                       
        for(int i=0; i<=N; i++){
            numericalFlux[i*(N+1)+N]    = 0.5*( boundaryRight[v][i*(N+1)]    + neighboringRight[v][i*(N+1)] ) ;   
            numericalFlux[i*(N+1)]    = 0.5*( boundaryLeft[v][i*(N+1)]     + neighboringLeft[v][i*(N+1)] ) ;  
        }
        /// vDash = -0.5*dy*D*v
        cblas_dgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dy, derivativeMatrix_x, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

        /// Adding the numeical Flux terms as necessary.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dy, fluxMatrix_right,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dy, fluxMatrix_left,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

        /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);

        delete[] numericalFlux;
        delete[] auxillaryVariable;
    }

    else if(fluxType == "rusanov") {
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        zeros(numericalFlux, (N+1)*(N+1)); 

        for(int i=0; i<=N; i++){
          // Normals nx, ny of the cell have been incorporated into the signs, need to set them separately !!
            numericalFlux[i*(N+1)+N] = 0.5*(boundaryRight[v][i*(N+1)] + neighboringRight[v][i*(N+1)] + MAX(fabs(boundaryRight[fluxVariable][i*(N+1)]), fabs(neighboringRight[fluxVariable][i*(N+1)]))*(boundaryRight[conserVar][i*(N+1)] - neighboringRight[conserVar][i*(N+1)])  ) ;   
            numericalFlux[i*(N+1)]   = 0.5*(boundaryLeft[v][i*(N+1)]  + neighboringLeft[v][i*(N+1)]  - MAX(fabs(boundaryLeft[fluxVariable][i*(N+1)]), fabs(neighboringLeft[fluxVariable][i*(N+1)]))*(boundaryLeft[conserVar][i*(N+1)] - neighboringLeft[conserVar][i*(N+1)])  ) ;   
        }

        /// vDash = -0.5*dy*D*v
        cblas_dgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dy, derivativeMatrix_x, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

        /// Adding the numeical Flux terms as necessary.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dy, fluxMatrix_right,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dy, fluxMatrix_left,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

        /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);

        delete[] numericalFlux;
        delete[] auxillaryVariable;

    }
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to get the variable `v` differentiated partially w.r.t. `y` and then store it in the
 * variable `vDash`. The function also takes `fluxType` as an input which would describe the numerical scheme that
 * should be used in order to obtain the derivative.
 *
 * @Param v         Variable which is to be differentiated.
 * @Param vDash     Variable in which the derivative is to be stored.
 * @Param conserVar Corresponding Conservative variable.
 * @Param fluxType  The type of flux that is to be used. eg "central"
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::delByDelY(int v, int vDash, int conserVar, string fluxType, int fluxVariable = 999) {
    double dy = (y_end - y_start);
    double dx = (x_end - x_start);

    if(fluxType == "central") {
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        zeros(numericalFlux, (N+1)*(N+1));                                                       
        for(int i=0; i<=N; i++){
            numericalFlux[N*(N+1)+i]    = 0.5*( boundaryTop[v][i]    + neighboringTop[v][i] ) ;   
            numericalFlux[i]            = 0.5*( boundaryBottom[v][i]     + neighboringBottom[v][i] ) ;  
        }
        /// vDash = -0.5*dy*D*v
        cblas_dgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dx, derivativeMatrix_y, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

        /// Adding the numeical Flux terms as necessary.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dx, fluxMatrix_top,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dx, fluxMatrix_bottom,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

        /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);

        delete[] numericalFlux;
        delete[] auxillaryVariable;
    }
    
    else if(fluxType == "rusanov") {
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        zeros(numericalFlux, (N+1)*(N+1));                                                       
        for(int i=0; i<=N; i++){
          // Normals nx, ny have been incorporated into the signs, need to set them separately !!
            numericalFlux[N*(N+1)+i]= 0.5*(boundaryTop[v][i] + neighboringTop[v][i] + MAX(fabs(boundaryTop[fluxVariable][i]), fabs(neighboringTop[fluxVariable][i]))*(boundaryTop[conserVar][i] - neighboringTop[conserVar][i]));
            numericalFlux[i]        = 0.5*(boundaryBottom[v][i] + neighboringBottom[v][i] - MAX(fabs(boundaryBottom[fluxVariable][i]), fabs(neighboringBottom[fluxVariable][i]))*(boundaryBottom[conserVar][i] - neighboringBottom[conserVar][i])); 
        }
        /// vDash = -0.5*dy*D*v
        cblas_dgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dx, derivativeMatrix_y, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

        /// Adding the numeical Flux terms as necessary.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dx, fluxMatrix_top,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dx, fluxMatrix_bottom,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

        /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);

        delete[] numericalFlux;
        delete[] auxillaryVariable;
    }
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to get the variable `v` differentiated partially w.r.t. `x` and then store it in the
 * variable `vDash`. The function also takes `fluxType` as an input which would describe the numerical scheme that
 * should be used in order to obtain the derivative.
 *
 * @Param v         Variable which is to be differentiated.
 * @Param vDash     Variable in which the derivative is to be stored.
 * @Param fluxType  The type of flux that is to be used. eg "central"
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::delByDelX(int v, int vDash, string fluxType) {
    double dy = (y_end - y_start);
    double dx = (x_end - x_start);
    
    if(fluxType == "central") {
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        zeros(numericalFlux, (N+1)*(N+1));                                                       
        for(int i=0; i<=N; i++){
            numericalFlux[i*(N+1)+N]    = 0.5*( boundaryRight[v][i*(N+1)]    + neighboringRight[v][i*(N+1)] ) ;   
            numericalFlux[i*(N+1)]    = 0.5*( boundaryLeft[v][i*(N+1)]     + neighboringLeft[v][i*(N+1)] ) ;  
        }
        /// vDash = -0.5*dy*D*v
        cblas_dgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dy, derivativeMatrix_x, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

        /// Adding the numeical Flux terms as necessary.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dy, fluxMatrix_right,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dy, fluxMatrix_left,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

        /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);

        delete[] numericalFlux;
        delete[] auxillaryVariable;
    }

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to get the variable `v` differentiated partially w.r.t. `y` and then store it in the
 * variable `vDash`. The function also takes `fluxType` as an input which would describe the numerical scheme that
 * should be used in order to obtain the derivative.
 *
 * @Param v         Variable which is to be differentiated.
 * @Param vDash     Variable in which the derivative is to be stored.
 * @Param fluxType  The type of flux that is to be used. eg "central"
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::delByDelY(int v, int vDash, string fluxType) {
    double dy = (y_end - y_start);
    double dx = (x_end - x_start);

    if(fluxType == "central") {
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        zeros(numericalFlux, (N+1)*(N+1));                                                       
        for(int i=0; i<=N; i++){
            numericalFlux[N*(N+1)+i]    = 0.5*( boundaryTop[v][i]    + neighboringTop[v][i] ) ;   
            numericalFlux[i]            = 0.5*( boundaryBottom[v][i]     + neighboringBottom[v][i] ) ;  
        }
        /// vDash = -0.5*dy*D*v
        cblas_dgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dx, derivativeMatrix_y, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

        /// Adding the numeical Flux terms as necessary.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dx, fluxMatrix_top,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dx, fluxMatrix_bottom,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

        /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
        cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);

        delete[] numericalFlux;
        delete[] auxillaryVariable;
    }
    
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function takes Mass Matrix as an input
 *
 * @Param m This is the massMatix array(actually a matrix, but implemented as a 1-d array).
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setMassMatrix(double *m) {
    massMatrix = m;
    return ;
}

void DG_Element_2d::setInverseMassMatrix(double* im) {
    inverseMassMatrix = im;
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function takes VanderMand Matrix as an input
 *
 * @Param vm This is the VanderMandMatrix array(actually a matrix, but implemented as a 1-d array).
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setVanderMandMatrix(double *vm) {
    vanderMandMatrix = vm;
    return ;
}

void DG_Element_2d::setInverseVanderMandMatrix(double* ivm) {
    inverseVanderMandMatrix = ivm;
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to set the x-derivative matrix.
 *
 * @Param d The array of the x-derivative matrix which is given as an input.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setderivateMatrix_x(double *d) {
    derivativeMatrix_x = d;
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to set the y-derivative matrix.
 *
 * @Param d The array of the y-derivative matrix which is given as an input.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setderivateMatrix_y(double *d) {
    derivativeMatrix_y = d;
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to set the transpose of x-derivative matrix.
 *
 * @Param d The array of transpose of x-derivative matrix which is given as an input.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setTransposederivateMatrix_x(double *id) {
    transposeDerivateMatrix_x = id;
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to set the transpose of y-derivative matrix.
 *
 * @Param d The array of transpose of y-derivative matrix which is given as an input.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setTransposederivateMatrix_y(double *id) {
    transposeDerivateMatrix_y = id;
    return ;
}

/***THE FOLLOWING 4 FUNCTIONS HAVE THE SIMILAR JOB TO SET THE FLUX MATRICES.****/

void DG_Element_2d::setTopFluxMatrix(double* f) {
    fluxMatrix_top = f;
    return ;
}
void DG_Element_2d::setRightFluxMatrix(double* f){
    fluxMatrix_right = f;
    return ;
}

void DG_Element_2d::setLeftFluxMatrix(double* f){
    fluxMatrix_left = f;
    return ;
}

void DG_Element_2d::setBottomFluxMatrix(double* f){
    fluxMatrix_bottom = f;
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @brief      This is used to extern the function `saxpy` to the class DG_Element_2d. $\mathbf{y} \mapsto a\mathbf{x} + \mathbf{y}$
 *
 * @param[in]  a     The coefficient `a` of `x`
 * @param[in]  x     The column vector `x`
 * @param[in]  y     The column vector `y`
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::axpy(double a, int x, int y) {
    cblas_daxpy((N+1)*(N+1), a, variable[x], 1, variable[y], 1);
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @brief      This is used to extern the function `sscal` to the class DG_Element_2d. $\mathbf{x} \mapsto a\mathbf{x}$
 *
 * @param[in]  a     The coefficient `a` of `x`
 * @param[in]  x     The column vector `x`
 * @param[in]  y     The column vector `y`
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::scal(double a, int x) {
    cblas_dscal((N+1)*(N+1), a, variable[x], 1);
    return ;
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
void DG_Element_2d::setFunctionsForVariables(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, unsigned, double*)> f, int z) {
    f(a, variable[x], b, variable[y], 1, (N+1)*(N+1), variable[z]);

    return ;
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
void DG_Element_2d::setFunctionsForVariables(double a, int w, double b, int x, double c, int y, function<void(double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z) {
    f(a, variable[w], b, variable[x], c, variable[y], 1, (N+1)*(N+1), variable[z]);

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
void DG_Element_2d::setFunctionsForVariables(double t, int a, double u, int b, double v, int c, double x, int d, function<void(double, double*, double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z) {
    f(t, variable[a], u, variable[b], v, variable[c], x, variable[d], 1, (N+1)*(N+1), variable[z]);

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
void DG_Element_2d::setFunctionsForBoundaryVariables(double a, int x, double b, int y, function<void(double, double*, double, double*, unsigned, unsigned, double*)> f, int z) {
    /*for(int i = 0 ; i <= N ; i++){
         boundaryTop[z][i] = f(boundaryTop[x][i], boundaryTop[y][i]);
         boundaryRight[z][i*(N+1)] = f(boundaryRight[x][i*(N+1)], boundaryRight[y][i*(N+1)]);
         boundaryLeft[z][i*(N+1)] = f(boundaryLeft[x][i*(N+1)], boundaryLeft[y][i*(N+1)]);
         boundaryBottom[z][i] = f(boundaryBottom[x][i], boundaryBottom[y][i]);
    }*/
    f(a, boundaryTop[x], b, boundaryTop[y], 1, N+1, boundaryTop[z]);
    f(a, boundaryBottom[x], b, boundaryBottom[y], 1, N+1, boundaryBottom[z]);
    f(a, boundaryRight[x], b, boundaryRight[y], N+1, N+1, boundaryRight[z]);
    f(a, boundaryLeft[x], b, boundaryLeft[y], N+1, N+1, boundaryLeft[z]);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of  Boundary variable z to f(w, x, y).
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
void DG_Element_2d::setFunctionsForBoundaryVariables(double a, int w, double b, int x, double c, int y, function<void(double, double*, double, double*, double, double*, unsigned, unsigned, double*)> f, int z) {
    /*for(int i = 0 ; i <= N ; i++){
         boundaryTop[z][i] = f(boundaryTop[w][i], boundaryTop[x][i], boundaryTop[y][i]);
         boundaryRight[z][i*(N+1)] = f(boundaryRight[w][i*(N+1)], boundaryRight[x][i*(N+1)], boundaryRight[y][i*(N+1)]);
         boundaryLeft[z][i*(N+1)] = f(boundaryLeft[w][i*(N+1)], boundaryLeft[x][i*(N+1)], boundaryLeft[y][i*(N+1)]);
         boundaryBottom[z][i] = f(boundaryBottom[w][i], boundaryBottom[x][i], boundaryBottom[y][i]);
    }*/
    f(a, boundaryTop[w], b, boundaryTop[x], c, boundaryTop[y], 1, N+1, boundaryTop[z]);
    f(a, boundaryBottom[w], b, boundaryBottom[x], c, boundaryBottom[y], 1, N+1, boundaryBottom[z]);
    f(a, boundaryRight[w], b, boundaryRight[x], c, boundaryRight[y], N+1, N+1, boundaryRight[z]);
    f(a, boundaryLeft[w], b, boundaryLeft[x], c, boundaryLeft[y], N+1, N+1, boundaryLeft[z]);

    return;
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
/*void DG_Element_2d::setFunctionsForVariablesCellCentered(int x, int y, function<double(double, double, double)> f, int z) {
    
    for(int i = 0 ; i < (N+1)*(N+1); i++)
        *variable[z] = f(variable[x][i],variable[y][i], *variable[z]);
    return ;
}*/



double DG_Element_2d::l2Norm(int v1, int v2) {
    double* diff = new double[(N+1)*(N+1)];
    cblas_dscal((N+1)*(N+1), 0.0, diff, 1);
    cblas_daxpy((N+1)*(N+1),  1.0, variable[v1], 1, diff, 1);
    cblas_daxpy((N+1)*(N+1), -1.0, variable[v2], 1, diff, 1);
    double norm2 =  (cblas_dnrm2((N+1)*(N+1),diff, 1));    
    delete[] diff;
    return norm2;
}
/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the maximum value in the Field .
 *
 * @Param v This is a int which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
/*double DG_Element_2d::FindMax(int v) {
    double Max = variable[v][0]; 
    for(int i = 0; i < (N+1)*(N+1); i++)
            Max = max(variable[v][i], Max);
    return Max;
}*/

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the minimum dx in the field.
 *
 * @Param dx To store the minimum dx in the element
 */
/* ----------------------------------------------------------------------------*/
/*void DG_Element_2d::FindMindx(int dx) {
    *variable[dx] = min(dxMin, dyMin);
    return ;
}*/

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the minimum dt in the field.
 *
 * @Param dt To find th mimimum dt from all dt.
 */
/* ----------------------------------------------------------------------------*/
/*double DG_Element_2d::FindMindt(int dt) {
   return *variable[dt];
}*/

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function finds the minimum dt in the element.
 *
 * @Param dx To store the minimum dx in the element
 * @Param dt Stores the minimum dt in an element
 * @Param U Stores the Maximum eigne value in the element
 * @Param CFL CFL to be used
 */
/* ----------------------------------------------------------------------------*/
/*void DG_Element_2d::FindTimestep(int dt, int dx, int U, double CFL) {
    *variable[dt] = (CFL/(0.0*N + 1.0)) * (*variable[dx])/(*variable[U]);
    return ;
}*/

/* ----------------------------------------------------------------------------*/
//--------------------------VIRTUAL FUNCTIONS----------------------------------//
/*-----------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates maps from a variable to its corresponding Boundary Setting methods
 *
 * @Param type This is the type of the boundary.
 * @Param b This defines the location of the boundary element.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::assignBoundary( string type, char b) {
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates Variables at the Boundary.
 *
 * @Param v This is the variable ,whose boundary values are to be fixed.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::updateBoundaryVariables(int v) {
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function performs MinMod for Boundary elements, eliminating all non-existing 
 * neighbouring cells.
 *
 * @Param m This int represents the moment of the variable
 * @Param Index This is represents to be limited.
 * @Param Alpha This is the scaling factor used in Lilia's Moment Limiter.
 * @Param R This is the pointer to RightNeighbour.
 * @Param L This is the pointer to LeftNeighbour.
 * @Param T This is the pointer to TopNeighbour.
 * @Param B This is the pointer to BottomNeighbour.
*/
/* ----------------------------------------------------------------------------*/
double DG_Element_2d::BoundaryMinMod(int m, int Index, double Alpha, DG_Element_2d* Right, DG_Element_2d* Left, DG_Element_2d* Top, DG_Element_2d* Bottom) {
    return 0.0 ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates Boundary values after a timestep.
 *
  * @Param time the solution time.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::updateBoundary(double time) {
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function adds conservative variables.
 *
 * @Param v This is the conservative variable 
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::addConservativeVariables(int v) {
    return ;
}

void DG_Element_2d::addConservativeVariables(vector<int> v) {
    return ;
}