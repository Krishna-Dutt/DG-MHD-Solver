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
  
  double *Var = variable[v];
  MaxVariable = Var[0];
  for (int i = 0; i < (N+1)*(N+1) ; ++i) {
    MaxVariable = MAX(MaxVariable,Var[i]);
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
  
  /*for(int b=0; b <(N+1)*(N+1); ++b) {
      variable[35][b] = Marker; // CellMarkerG
  }*/
 
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
    double *VarM = variable[m], *VarModM = variable[modm]; 
    double *RNVarM = rightNeighbor->variable[m];
    double *LNVarM = leftNeighbor->variable[m];
    double *TNVarM = topNeighbor->variable[m];
    double *BNVarM = bottomNeighbor->variable[m];
    
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

       // Original minmod detector
       Temp1 = MinMod(VarM[Tempi], AlphaN*(RNVarM[Tempi-1] -VarM[Tempi-1]), AlphaN*(VarM[Tempi-1] -LNVarM[Tempi-1]) , AlphaN*(TNVarM[Tempi-(N+1)] -VarM[Tempi-(N+1)]), AlphaN*(VarM[Tempi-(N+1)] -BNVarM[Tempi-(N+1)]));
       Temp2 = MinMod(VarM[Tempj], AlphaN*(RNVarM[Tempj-1] -VarM[Tempj-1]), AlphaN*(VarM[Tempj-1] -LNVarM[Tempj-1]) , AlphaN*(TNVarM[Tempj-(N+1)] -VarM[Tempj-(N+1)]), AlphaN*(VarM[Tempj-(N+1)] -BNVarM[Tempj-(N+1)]));
       
       if (abs(Temp1-VarModM[Tempi]) > epsilon || abs(Temp2-VarModM[Tempj]) > epsilon ) {
         VarModM[Tempi] = Temp1;
         VarModM[Tempj] = Temp2;
       }
       else 
       {
         return ; // Need to exit both loops
       }
     }
     // Special Case for end values, when Tempi or Tempj access zero order polynomials !!
       Tempi = i-j;
       Tempj = i - j*(N+1);

       // Original Detector
       Temp1 = MinMod(VarM[Tempi], AlphaN*(TNVarM[Tempi-(N+1)] -VarM[Tempi-(N+1)]), AlphaN*(VarM[Tempi-(N+1)] -BNVarM[Tempi-(N+1)]));
       Temp2 = MinMod(VarM[Tempj], AlphaN*(RNVarM[Tempj-1] -VarM[Tempj-1]), AlphaN*(VarM[Tempj-1] -LNVarM[Tempj-1]));
       
       if ( abs(Temp1-VarModM[Tempi]) > epsilon || abs(Temp2-VarModM[Tempj]) > epsilon ) {
         VarModM[Tempi] = Temp1;
         VarModM[Tempj] = Temp2;
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
 * @Synopsis  This function computes the moments for a given variable using VandeMandMatrix with Legendre Basis.
 *
 * @Param v This is the variable whose moments are to be captured.
 * @Param m This is the variable used to store the computed moments.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::computeMoments(int *v, int *m, unsigned size) {
  /// Multiplying inverse of VanderMand Matrix with the variable array to get the corresponding moments.
  for(int i=0; i< size; ++i) {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 1.0, inverseVanderMandMatrix,(N+1)*(N+1), variable[v[i]],1,0,variable[m[i]],1);
  }
  
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
void DG_Element_2d::convertMomentToVariable(int *m, int *v, unsigned size) {
  /// Multiplying  VanderMand Matrix with the moments to obtained the nodal values of the variable.
  for(int i=0; i< size; ++i) {
      cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 1.0, vanderMandMatrix,(N+1)*(N+1), variable[m[i]],1,0,variable[v[i]],1);
  }
  
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
void DG_Element_2d::limitMoments(int *M, int *Modm, unsigned Index, unsigned size ) {

   // Checking if cell marker is not equal to zero
    int count, Tempi, Tempj, i, j, m, modm;
    count = N+1;
    double Temp1, Temp2, AlphaN;
    double epsilon = 1e-13;
    AlphaN = sqrt((2.0*N -1.0)/(2.0*N +1));  // Similar to a diffusion coefficient

    if ( N == 1) {
        epsilon = 0.0 ;
    }
    else {
        epsilon = 1e-16;
    }
    int counter = 0;
    for(int temp = 0 ; temp<size; ++temp) {
        m = M[temp]; modm = Modm[temp];
        double *VarM = variable[m], *VarModM = variable[modm]; 
        double *RNVarM = rightNeighbor->variable[m];
        double *LNVarM = leftNeighbor->variable[m];
        double *TNVarM = topNeighbor->variable[m];
        double *BNVarM = bottomNeighbor->variable[m];
    
        // Ensuring that Cell avergae remains the  same after limiting !!
        variable[modm][0] = variable[m][0];

        count = N+1;
        counter = 0;
        AlphaN = sqrt((2.0*N -1.0)/(2.0*N +1));
        for(i=Index; i > 0 && counter == 0; i = i - (N+2)) {
          --count;
          //AlphaN = sqrt((2.0*(count)-1.0)/(2.0*(count)+1.0)); 
          //AlphaN = 0.5*sqrt((4.0*(count)-1.0)/(2.0*(count)+1.0));
          //AlphaN = 0.5/sqrt(4.0*count*count -1.0);
          for(j=0; j < count && counter == 0; ++j) {
             Tempi = i-j;
             Tempj = i - j*(N+1);
             // Original minmod detector
             Temp1 = MinMod(VarM[Tempi], AlphaN*(RNVarM[Tempi-1] -VarM[Tempi-1]), AlphaN*(VarM[Tempi-1] -LNVarM[Tempi-1]) , AlphaN*(TNVarM[Tempi-(N+1)] -VarM[Tempi-(N+1)]), AlphaN*(VarM[Tempi-(N+1)] -BNVarM[Tempi-(N+1)]));
             Temp2 = MinMod(VarM[Tempj], AlphaN*(RNVarM[Tempj-1] -VarM[Tempj-1]), AlphaN*(VarM[Tempj-1] -LNVarM[Tempj-1]) , AlphaN*(TNVarM[Tempj-(N+1)] -VarM[Tempj-(N+1)]), AlphaN*(VarM[Tempj-(N+1)] -BNVarM[Tempj-(N+1)]));
       
             if (abs(Temp1-VarModM[Tempi]) > epsilon || abs(Temp2-VarModM[Tempj]) > epsilon ) {
                VarModM[Tempi] = Temp1;
                VarModM[Tempj] = Temp2;
            }
             else 
             {
               counter = 1; // Need to exit both loops
             }
          }
          // Special Case for end values, when Tempi or Tempj access zero order polynomials !!
          if(counter == 0) {
                Tempi = i-j;
                Tempj = i - j*(N+1);
                // Original Detector
                Temp1 = MinMod(VarM[Tempi], AlphaN*(TNVarM[Tempi-(N+1)] -VarM[Tempi-(N+1)]), AlphaN*(VarM[Tempi-(N+1)] -BNVarM[Tempi-(N+1)]));
                Temp2 = MinMod(VarM[Tempj], AlphaN*(RNVarM[Tempj-1] -VarM[Tempj-1]), AlphaN*(VarM[Tempj-1] -LNVarM[Tempj-1]));
       
                if ( abs(Temp1-VarModM[Tempi]) > epsilon || abs(Temp2-VarModM[Tempj]) > epsilon ) {
                  VarModM[Tempi] = Temp1;
                  VarModM[Tempj] = Temp2;
                }
                else //if( Temp1 !=0 && Temp2 !=0)
                 {
                   counter = 1; // Need to exit both loops
                 }
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
 * @Param array C moment of characteristic variables, to store modified moments.
 * @Param Index Index correspoding to Characteristic variable.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::limitMoments(int *V, int *C, unsigned Index) {
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
    
    int counter = 0;
    for(int temp=0; temp< Dimension; ++temp) {
        count = N+1;
        counter = 0;
        double *VarC = variable[C[temp]];

        for(i=(N+1)*(N+1)-1; (i > 0) && (counter == 0); i = i - (N+2)) {
          --count;

          AlphaN = sqrt((2.0*(count)-1.0)/(2.0*(count)+1.0));
          //AlphaN = 0.5/sqrt(4.0*count*count -1.0);
          //AlphaN = 0.25*(4.0*count-1.0)/sqrt(4.0*count*count -1.0);
          for(j=0; (j < count) && (counter == 0); ++j) {
             Tempi = i-j;
             Tempj = i - j*(N+1);

             Sum1 = Sum2 = Sum3 = Sum4 = 0;
             for(int z=0; z<Dimension; ++z) {
                 Sum1 += LeftEigenMatrix[temp*Dimension+z] * rightNeighbor->variable[V[z]][Tempi-1];
                 Sum2 += LeftEigenMatrix[temp*Dimension+z] * leftNeighbor->variable[V[z]][Tempi-1]; 
                 Sum3 += LeftEigenMatrix[temp*Dimension+z] * topNeighbor->variable[V[z]][Tempi-(N+1)];
                 Sum4 += LeftEigenMatrix[temp*Dimension+z] * bottomNeighbor->variable[V[z]][Tempi-(N+1)];   
             }

             Temp1 = MinMod(VarC[Tempi], AlphaN*(Sum1 -VarC[Tempi-1]), AlphaN*(VarC[Tempi-1] -Sum2), AlphaN*(Sum3 -VarC[Tempi-(N+1)]), AlphaN*(VarC[Tempi-(N+1)] -Sum4));
                          
             Sum1 = Sum2 = Sum3 = Sum4 = 0;
             for(int z=0; z<Dimension; ++z) {
                 Sum1 += LeftEigenMatrix[temp*Dimension+z] * rightNeighbor->variable[V[z]][Tempj-1];
                 Sum2 += LeftEigenMatrix[temp*Dimension+z] * leftNeighbor->variable[V[z]][Tempj-1]; 
                 Sum3 += LeftEigenMatrix[temp*Dimension+z] * topNeighbor->variable[V[z]][Tempj-(N+1)];
                 Sum4 += LeftEigenMatrix[temp*Dimension+z] * bottomNeighbor->variable[V[z]][Tempj-(N+1)]; 
             }

             Temp2 = MinMod(VarC[Tempj], AlphaN*(Sum1 -VarC[Tempj-1]), AlphaN*(VarC[Tempj-1] -Sum2), AlphaN*(Sum3 -VarC[Tempj-(N+1)]), AlphaN*(VarC[Tempj-(N+1)] -Sum4));

             if (abs(Temp1-VarC[Tempi]) > epsilon || abs(Temp2-VarC[Tempj]) > epsilon ) {
                 VarC[Tempi] = Temp1;
                 VarC[Tempj] = Temp2;
             }
             else {
                 counter = 1.0 ;
             }
          } 
             // Special Case for end values, when Tempi or Tempj access zero order polynomials !!
           if(counter == 0) {
             Tempi = i-j;
             Tempj = i - j*(N+1);
    
             Sum1 = Sum2 = Sum3 = Sum4 = 0;
             for(int z=0; z<Dimension; ++z) {
                 Sum1 += LeftEigenMatrix[temp*Dimension+z] * topNeighbor->variable[V[z]][Tempi-(N+1)]; 
                 Sum2 += LeftEigenMatrix[temp*Dimension+z] * bottomNeighbor->variable[V[z]][Tempi-(N+1)]; 
                 Sum3 += LeftEigenMatrix[temp*Dimension+z] * rightNeighbor->variable[V[z]][Tempj-1]; 
                 Sum4 += LeftEigenMatrix[temp*Dimension+z] * leftNeighbor->variable[V[z]][Tempj-1]; 
             }
            
             Temp1 = MinMod(VarC[Tempi], AlphaN*(Sum1 -VarC[Tempi-(N+1)]), AlphaN*(VarC[Tempi-(N+1)] -Sum2));
             Temp2 = MinMod(VarC[Tempj], AlphaN*(Sum3 -VarC[Tempj-1]), AlphaN*(VarC[Tempj-1] -Sum4));
             
             if ( abs(Temp1-VarC[Tempi]) > epsilon || abs(Temp2-VarC[Tempj]) > epsilon ) {
                VarC[Tempi] = Temp1;
                VarC[Tempj] = Temp2;
             }
             else {
                 counter = 1.0 ;
             }
           }  
             
        }

    }

  return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function which sets the Right and Left Eigen Matrix corresponding to directionalong flow
 * 
 * @Param dimension No of system of Equations
 * @Param REMp Pointer to Right Eigen Matrix
 * @Param LEMp Pointer to Left Eigen Matrix
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setEigenMatrices(unsigned dimension, double *REMp, double *LEMp) {
  
 RightEigenMatrix = REMp;
 LeftEigenMatrix = LEMp;

 Dimension = dimension;

  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to find  the Right and Left Eigen Matrix corresponding to direction along flow
 * for Euler system
 * 
 * @Param array V array of index corresponding to field variables required to compute eigen matrices.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::findEigenMatricesEuler(int *V) {
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
  
 
  return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to find  the Right and Left Eigen Matrix corresponding to direction along flow
 * for MHD system
 * 
 * @Param array V array of index corresponding to field variables required to compute eigen matrices.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::findEigenMatricesMHD(int *V) {
  double u, v, w, c, H, P;
  double epsilon = 1e-10;
  double q, qn, Bn, nx, ny, qe, cfn, can, csn, Pt, BdotB, UdotB, ch, Bcross, Ucross;
  double DVx, DVy, DVz, D, DE, Bx, By, Bz, Si, dPdx, dPdy;
  double theta1, theta2, scale;
  D   = 0.5 * variable[V[0]][0];
  DVx = 0.5 * variable[V[1]][0];
  DVy = 0.5 * variable[V[2]][0];
  DVz = 0.5 * variable[V[3]][0];
  DE  = 0.5 * variable[V[4]][0];
  Bx = 0.5 * variable[V[5]][0];
  By = 0.5 * variable[V[6]][0];
  Bz = 0.5 * variable[V[7]][0];
  Si = 0.5 * variable[V[8]][0];
  dPdx = 0.5 * variable[V[4]][0];
  dPdy = 0.5 * variable[V[5]][0];
  
  u = DVx/D;
  v = DVy/D;
  w = DVz/D;
  
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
  Bn = Bx*nx + By*ny;
  BdotB = Bx*Bx + By*By + Bz*Bz;
  UdotB = u*Bx + v*By + w*Bz;
  Bcross = Bx*ny - By*nx;
  Ucross = u*ny - v*nx;
  q = sqrt(u*u + v*v + w*w);
  P = (gamma-1.0) * ( DE - 0.5*q*q*D - 0.5*BdotB ); 
  Pt = P + 0.5 * BdotB;
  
  c = (gamma*P/D);
  can = abs(Bn)/D;
  cfn = sqrt( 0.5 * ( c + BdotB/D  + sqrt( pow(c + BdotB/D, 2.0) - 4.0*c*Bn*Bn/D ) ) );
  csn = sqrt( 0.5 * ( c + BdotB/D  - sqrt( pow(c + BdotB/D, 2.0) - 4.0*c*Bn*Bn/D ) ) );

  
  // Checks for Zero Magnetic field !!
  if ( (abs(Bx) + abs(By) + abs(Bz)) < epsilon ) 
  {
      zeros(RightEigenMatrix, Dimension);
      zeros(LeftEigenMatrix, Dimension);
  }
  // Setting Right Eigen Matrix
  // Lambda = qn 
  RightEigenMatrix[4] = 1.0; RightEigenMatrix[13] = u; RightEigenMatrix[22] = v; RightEigenMatrix[31] = w;
  RightEigenMatrix[40] = 0.5*q*q; RightEigenMatrix[49] =  RightEigenMatrix[58] = RightEigenMatrix[67] = RightEigenMatrix[76] = 0.0 ;
  
  // Condition for Bz and Bn
  // Lambda = qn + can
  RightEigenMatrix[6] = 0.0; RightEigenMatrix[15] = -Bz*ny; RightEigenMatrix[24] = Bz*nx;
  RightEigenMatrix[33] = Bcross; RightEigenMatrix[42] = -Bz*Ucross -w*Bcross;
  RightEigenMatrix[51] = can*ny*Bz/Bn; RightEigenMatrix[60] = -can*nx*Bz/Bn;
  RightEigenMatrix[69] = -can*Bcross/Bn; RightEigenMatrix[76] = 0.0;  

  // Lambda = qn - can
  RightEigenMatrix[2] = 0.0; RightEigenMatrix[11] = -Bz*ny; RightEigenMatrix[20] = Bz*nx;
  RightEigenMatrix[29] = Bcross; RightEigenMatrix[38] = -Bz*Ucross -w*Bcross;
  RightEigenMatrix[47] = -can*ny*Bz/Bn; RightEigenMatrix[56] = can*nx*Bz/Bn;
  RightEigenMatrix[65] = can*Bcross/Bn; RightEigenMatrix[74] = 0.0; 

  // Lambda = qn + csn
  theta2 = -Bcross*( (gamma-1.0)*Bn*Bn/D + (2.0-gamma)*csn*csn );
  theta1 = (1.0-gamma)*Bn * ( (Bn*UdotB/D + csn*(DE + Pt)/D)*nx*(pow(csn,2.0) - pow(Bn,2.0)/D) -nx*Bn*w*Bz*csn*csn/D -(Bn/D)*( (Bx*(qn + csn) -w*Bz*nx)*(pow(csn,2.0)-pow(Bn,2.0)/D) - nx*csn*Bz*Bz*Bn/D) ) 
           - csn* ( (gamma-1.0)*(2.0*w*w - 0.5*q*q)*Bn*nx*(pow(csn,2.0) - pow(Bn,2.0)/D) + (gamma-1.0)*(nx*w*Bn*Bn*Bz*csn/D) + (2.0-gamma)*( -Bcross*(csn*csn*ny + csn*v)*(pow(csn,2.0)-pow(Bn,2.0)/D) + nx*Bz*Bn*csn*csn*Bz/D ) );
  RightEigenMatrix[5] = theta2*nx*(pow(csn,2.0) - pow(Bn,2.0)/D); RightEigenMatrix[14] = (theta2*(qn + csn) + theta1*ny)*(pow(csn,2.0) - pow(Bn,2.0)/D);
  RightEigenMatrix[23] = -theta1*nx*(pow(csn,2.0) - pow(Bn,2.0)/D); RightEigenMatrix[32] = -theta2*nx*(w*(pow(csn,2.0) - pow(Bn,2.0)/D) + csn*Bz*Bn/D);
  RightEigenMatrix[59] = nx*(pow(csn,2.0) - pow(Bn,2.0)/D)*((csn*csn*ny + csn*v)*theta2 + csn*theta1)/Bn; RightEigenMatrix[50] = -ny*(pow(csn,2.0) - pow(Bn,2.0)/D)*((csn*csn*ny + csn*v)*theta2 + csn*theta1)/Bn; 
  RightEigenMatrix[68] = nx*theta2*csn*csn*Bz/D; RightEigenMatrix[77] = 0.0;
  RightEigenMatrix[41] = ((0.5*(gamma-1.0)*q*q - csn*csn)*RightEigenMatrix[5] + (1-gamma)*(u*RightEigenMatrix[14] + v*RightEigenMatrix[23] + w*RightEigenMatrix[32]) + (2.0-gamma)*( Bx*RightEigenMatrix[50] + By*RightEigenMatrix[59] + Bz*RightEigenMatrix[68]))/(1.0-gamma) ;
  
  // Lambda = qn - csn
  theta2 = -Bcross*( (gamma-1.0)*Bn*Bn/D + (2.0-gamma)*csn*csn );
  theta1 = (1.0-gamma)*Bn * ( (Bn*UdotB/D - csn*(DE + Pt)/D)*nx*(pow(csn,2.0) - pow(Bn,2.0)/D) -nx*Bn*w*Bz*csn*csn/D -(Bn/D)*( (Bx*(qn - csn) -w*Bz*nx)*(pow(csn,2.0)-pow(Bn,2.0)/D) + nx*csn*Bz*Bz*Bn/D) ) 
           + csn* ( (gamma-1.0)*(2.0*w*w - 0.5*q*q)*Bn*nx*(pow(csn,2.0) - pow(Bn,2.0)/D) + (gamma-1.0)*(-nx*w*Bn*Bn*Bz*csn/D) + (2.0-gamma)*( -Bcross*(csn*csn*ny - csn*v)*(pow(csn,2.0)-pow(Bn,2.0)/D) + nx*Bz*Bn*csn*csn*Bz/D ) );
  RightEigenMatrix[3] = theta2*nx*(pow(csn,2.0) - pow(Bn,2.0)/D); RightEigenMatrix[12] = (theta2*(qn - csn) + theta1*ny)*(pow(csn,2.0) - pow(Bn,2.0)/D);
  RightEigenMatrix[21] = -theta1*nx*(pow(csn,2.0) - pow(Bn,2.0)/D); RightEigenMatrix[30] = -theta2*nx*(w*(pow(csn,2.0) - pow(Bn,2.0)/D) - csn*Bz*Bn/D);
  RightEigenMatrix[57] = nx*(pow(csn,2.0) - pow(Bn,2.0)/D)*((csn*csn*ny - csn*v)*theta2 - csn*theta1)/Bn; RightEigenMatrix[48] = -ny*(pow(csn,2.0) - pow(Bn,2.0)/D)*((csn*csn*ny - csn*v)*theta2 - csn*theta1)/Bn; 
  RightEigenMatrix[66] = nx*theta2*csn*csn*Bz/D; RightEigenMatrix[75] = 0.0;
  RightEigenMatrix[39] = ((0.5*(gamma-1.0)*q*q - csn*csn)*RightEigenMatrix[3] + (1-gamma)*(u*RightEigenMatrix[12] + v*RightEigenMatrix[21] + w*RightEigenMatrix[30]) + (2.0-gamma)*( Bx*RightEigenMatrix[48] + By*RightEigenMatrix[57] + Bz*RightEigenMatrix[66]))/(1.0-gamma) ;
  
  // Lambda = qn + cfn
  theta2 = -Bcross*( (gamma-1.0)*Bn*Bn/D + (2.0-gamma)*cfn*cfn );
  theta1 = (1.0-gamma)*Bn * ( (Bn*UdotB/D + cfn*(DE + Pt)/D)*nx*(pow(cfn,2.0) - pow(Bn,2.0)/D) -nx*Bn*w*Bz*cfn*cfn/D -(Bn/D)*( (Bx*(qn + cfn) -w*Bz*nx)*(pow(cfn,2.0)-pow(Bn,2.0)/D) - nx*cfn*Bz*Bz*Bn/D) ) 
           - cfn* ( (gamma-1.0)*(2.0*w*w - 0.5*q*q)*Bn*nx*(pow(cfn,2.0) - pow(Bn,2.0)/D) + (gamma-1.0)*(nx*w*Bn*Bn*Bz*cfn/D) + (2.0-gamma)*( -Bcross*(cfn*cfn*ny + cfn*v)*(pow(cfn,2.0)-pow(Bn,2.0)/D) + nx*Bz*Bn*cfn*cfn*Bz/D ) );
  RightEigenMatrix[7] = theta2*nx*(pow(cfn,2.0) - pow(Bn,2.0)/D); RightEigenMatrix[16] = (theta2*(qn + cfn) + theta1*ny)*(pow(cfn,2.0) - pow(Bn,2.0)/D);
  RightEigenMatrix[25] = -theta1*nx*(pow(cfn,2.0) - pow(Bn,2.0)/D); RightEigenMatrix[34] = -theta2*nx*(w*(pow(cfn,2.0) - pow(Bn,2.0)/D) + cfn*Bz*Bn/D);
  RightEigenMatrix[61] = nx*(pow(cfn,2.0) - pow(Bn,2.0)/D)*((cfn*cfn*ny + cfn*v)*theta2 + cfn*theta1)/Bn; RightEigenMatrix[52] = -ny*(pow(cfn,2.0) - pow(Bn,2.0)/D)*((cfn*cfn*ny + cfn*v)*theta2 + cfn*theta1)/Bn; 
  RightEigenMatrix[70] = nx*theta2*cfn*cfn*Bz/D; RightEigenMatrix[79] = 0.0;
  RightEigenMatrix[43] = ((0.5*(gamma-1.0)*q*q - cfn*cfn)*RightEigenMatrix[7] + (1-gamma)*(u*RightEigenMatrix[16] + v*RightEigenMatrix[25] + w*RightEigenMatrix[34]) + (2.0-gamma)*( Bx*RightEigenMatrix[52] + By*RightEigenMatrix[61] + Bz*RightEigenMatrix[70]))/(1.0-gamma) ;
  
  // Lambda = qn - cfn
  theta2 = -Bcross*( (gamma-1.0)*Bn*Bn/D + (2.0-gamma)*cfn*cfn );
  theta1 = (1.0-gamma)*Bn * ( (Bn*UdotB/D - cfn*(DE + Pt)/D)*nx*(pow(cfn,2.0) - pow(Bn,2.0)/D) -nx*Bn*w*Bz*cfn*cfn/D -(Bn/D)*( (Bx*(qn - cfn) -w*Bz*nx)*(pow(cfn,2.0)-pow(Bn,2.0)/D) + nx*cfn*Bz*Bz*Bn/D) ) 
           + cfn* ( (gamma-1.0)*(2.0*w*w - 0.5*q*q)*Bn*nx*(pow(cfn,2.0) - pow(Bn,2.0)/D) + (gamma-1.0)*(-nx*w*Bn*Bn*Bz*cfn/D) + (2.0-gamma)*( -Bcross*(cfn*cfn*ny - cfn*v)*(pow(cfn,2.0)-pow(Bn,2.0)/D) + nx*Bz*Bn*cfn*cfn*Bz/D ) );
  RightEigenMatrix[1] = theta2*nx*(pow(cfn,2.0) - pow(Bn,2.0)/D); RightEigenMatrix[10] = (theta2*(qn - cfn) + theta1*ny)*(pow(cfn,2.0) - pow(Bn,2.0)/D);
  RightEigenMatrix[19] = -theta1*nx*(pow(cfn,2.0) - pow(Bn,2.0)/D); RightEigenMatrix[28] = -theta2*nx*(w*(pow(cfn,2.0) - pow(Bn,2.0)/D) - cfn*Bz*Bn/D);
  RightEigenMatrix[55] = nx*(pow(cfn,2.0) - pow(Bn,2.0)/D)*((cfn*cfn*ny - cfn*v)*theta2 - cfn*theta1)/Bn; RightEigenMatrix[46] = -ny*(pow(cfn,2.0) - pow(Bn,2.0)/D)*((cfn*cfn*ny - cfn*v)*theta2 - cfn*theta1)/Bn; 
  RightEigenMatrix[64] = nx*theta2*cfn*cfn*Bz/D; RightEigenMatrix[73] = 0.0;
  RightEigenMatrix[37] = ((0.5*(gamma-1.0)*q*q - cfn*cfn)*RightEigenMatrix[1] + (1-gamma)*(u*RightEigenMatrix[10] + v*RightEigenMatrix[19] + w*RightEigenMatrix[28]) + (2.0-gamma)*( Bx*RightEigenMatrix[46] + By*RightEigenMatrix[55] + Bz*RightEigenMatrix[64]))/(1.0-gamma) ;
  
  // Lambda = ch
  theta1 = (pow(ch - qn,2.0)*(Bn*Bn - BdotB) + BdotB*( Bn*Bn/D - (ch - qn)))/D;
  theta2 = Bn*(ch - qn)*( Bn*Bn/D - BdotB/D + gamma*(pow(ch - qn,2.0) - Bn*Bn/D) ) + (qn*Bn - UdotB)*( -2.0*pow(ch - qn,2.0) + (1.0-gamma)*Bn*Bn/D );
  RightEigenMatrix[80] = ch * (pow(ch - qn, 2.0) - Bn*Bn/D)*theta1;
  RightEigenMatrix[71] = theta1*(Bz*Bn/D - w*(ch - qn)) - theta2*(pow(ch - qn, 2.0)*Bz/D); 
  RightEigenMatrix[62] = theta1*( -Bcross*nx*Bn/D + (ch - qn)*Ucross*nx + ny*( pow(ch - qn,2.0) - Bn*Bn/D )) - theta2*(pow(ch - qn,2.0)*(By - Bn*ny)/D);
  RightEigenMatrix[53] = theta1*( Bcross*ny*Bn/D - (ch - qn)*Ucross*ny + nx*( pow(ch - qn,2.0) - Bn*Bn/D )) - theta2*(pow(ch - qn,2.0)*Bcross*ny/D);
  RightEigenMatrix[8] = -theta2*(pow(ch - qn,2.0) - Bn*Bn/D) ;
  RightEigenMatrix[17] = theta1*( -(ch - qn)*ny*Bcross + Bn*ny*Ucross ) -theta2*( u*(pow(ch - qn,2.0) - Bn*Bn/D) + (ch - qn)*(pow(ch - qn,2.0)*nx - Bx*Bn/D) ) ;
  RightEigenMatrix[26] = theta1*( (ch - qn)*nx*Bcross - Bn*nx*Ucross ) -theta2*( v*(pow(ch - qn,2.0) - Bn*Bn/D) + (ch - qn)*(pow(ch - qn,2.0)*ny - By*Bn/D) ); 
  RightEigenMatrix[35] = theta1*( w*Bn - Bz*( ch - qn) ) - theta2*( w*(pow(ch - qn,2.0) - Bn*Bn/D) - Bz*Bn*( ch - qn)/D );
  RightEigenMatrix[44] = (( pow(ch - qn,2.0) - 0.5*(gamma-1.0)*q*q )*RightEigenMatrix[8] + (gamma-1.0)*(u*RightEigenMatrix[17] + v*RightEigenMatrix[26] + w*RightEigenMatrix[35])
                         + (gamma-2.0)*( Bx*RightEigenMatrix[53] + By*RightEigenMatrix[62] + Bz*RightEigenMatrix[71]) + (2*Bn/ch)*RightEigenMatrix[80] )/(gamma-1.0);  

  // Lambda = -ch
  theta1 = (pow(-ch - qn,2.0)*(Bn*Bn - BdotB) + BdotB*( Bn*Bn/D - (-ch - qn)))/D;
  theta2 = Bn*(-ch - qn)*( Bn*Bn/D - BdotB/D + gamma*(pow(-ch - qn,2.0) - Bn*Bn/D) ) + (qn*Bn - UdotB)*( -2.0*pow(-ch - qn,2.0) + (1.0-gamma)*Bn*Bn/D );
  RightEigenMatrix[72] = -ch * (pow(-ch - qn, 2.0) - Bn*Bn/D)*theta1;
  RightEigenMatrix[63] = theta1*(Bz*Bn/D - w*(-ch - qn)) - theta2*(pow(-ch - qn, 2.0)*Bz/D); 
  RightEigenMatrix[54] = theta1*( -Bcross*nx*Bn/D + (-ch - qn)*Ucross*nx + ny*( pow(-ch - qn,2.0) - Bn*Bn/D )) - theta2*(pow(-ch - qn,2.0)*(By - Bn*ny)/D);
  RightEigenMatrix[45] = theta1*( Bcross*ny*Bn/D - (-ch - qn)*Ucross*ny + nx*( pow(-ch - qn,2.0) - Bn*Bn/D )) - theta2*(pow(-ch - qn,2.0)*Bcross*ny/D);
  RightEigenMatrix[0] = -theta2*(pow(-ch - qn,2.0) - Bn*Bn/D) ;
  RightEigenMatrix[9] = theta1*( -(-ch - qn)*ny*Bcross + Bn*ny*Ucross ) -theta2*( u*(pow(-ch - qn,2.0) - Bn*Bn/D) + (-ch - qn)*(pow(-ch - qn,2.0)*nx - Bx*Bn/D) ) ;
  RightEigenMatrix[18] = theta1*( (-ch - qn)*nx*Bcross - Bn*nx*Ucross ) -theta2*( v*(pow(-ch - qn,2.0) - Bn*Bn/D) + (-ch - qn)*(pow(-ch - qn,2.0)*ny - By*Bn/D) ); 
  RightEigenMatrix[27] = theta1*( w*Bn - Bz*(-ch - qn) ) - theta2*( w*(pow(-ch - qn,2.0) - Bn*Bn/D) - Bz*Bn*(-ch - qn)/D );
  RightEigenMatrix[36] = (( pow(-ch - qn,2.0) - 0.5*(gamma-1.0)*q*q )*RightEigenMatrix[0] + (gamma-1.0)*(u*RightEigenMatrix[9] + v*RightEigenMatrix[18] + w*RightEigenMatrix[27])
                         + (gamma-2.0)*( Bx*RightEigenMatrix[45] + By*RightEigenMatrix[54] + Bz*RightEigenMatrix[63]) + (2*Bn/ch)*RightEigenMatrix[72] )/(gamma-1.0);  
  


  // Setting Left Eigen Matrix;
  //inverse(RightEigenMatrix,LeftEigenMatrix,Dimension*Dimension);
  // Lambda = qn
  if ( abs(Bcross) > epsilon )
   {
      LeftEigenMatrix[36] = ( -(DE + Pt)/D + BdotB/D + qn*qn ) * Bcross - Ucross * Bz * w;
      LeftEigenMatrix[37] = -(Bz*w*nx + qn *(Bx - Bn*nx)); LeftEigenMatrix[38] = (Bz*w*ny + qn*(By - Bn*ny));
      LeftEigenMatrix[39] = -Bcross*w; LeftEigenMatrix[40] = Bcross;
      LeftEigenMatrix[41] = -Bcross*Bcross*ny; LeftEigenMatrix[42] = Bcross*Bcross*nx;
      LeftEigenMatrix[43] = -Bz*Bcross; LeftEigenMatrix[44] = 0.0;

      scale = Bcross*( BdotB/D + 0.5*q*q - ( w*w + (DE + Pt)/D ) );
      for(int i=36; i<=44; ++i) LeftEigenMatrix[i] = LeftEigenMatrix[i]/scale;

  }

  else {
      LeftEigenMatrix[36] = q*q + (Bn*Bn + Bz*Bz - DE - Pt)/D ;
      LeftEigenMatrix[37] = -u; LeftEigenMatrix[38] = -v;
      LeftEigenMatrix[39] = -w; LeftEigenMatrix[40] = 1;
      LeftEigenMatrix[41] = 0; LeftEigenMatrix[42] = 0;
      LeftEigenMatrix[43] = -Bz; LeftEigenMatrix[44] = 0.0;

      scale = 0.5*q*q + (Bn*Bn + Bz*Bz - DE - Pt)/D;
      for(int i=36; i<=44; ++i) LeftEigenMatrix[i] = LeftEigenMatrix[i]/scale;

  }
  
  
  // Condition for Bn = 0.0
  // Lambda  = qn + can
  scale = 2.0*( Bn*(Bn*Bn - BdotB)/D + can*(UdotB - qn*Bn) );
  if ( abs(scale) < epsilon) scale = 1.0;
  theta1 = -Bn*Bcross/D + can*Ucross;
  LeftEigenMatrix[54] =  -( (-Ucross +can*Bcross/Bn)*(w*can -Bz*Bn/D) + (w - can*Bz/Bn)*theta1) ;  
  LeftEigenMatrix[55] = -ny*( w*can -Bz*Bn/D); LeftEigenMatrix[56] = nx*( w*can -Bz*Bn/D);
  LeftEigenMatrix[57] = theta1; LeftEigenMatrix[58] = 0.0;
  LeftEigenMatrix[59] = -ny*(-can*D/Bn)*(w*can - Bz*Bn/D); LeftEigenMatrix[60] = nx*(-can*D/Bn)*(w*can - Bz*Bn/D);
  LeftEigenMatrix[61] = theta1*(-can*D/Bn); LeftEigenMatrix[62] = 0.0;
  for(int i=54; i<=62; ++i) LeftEigenMatrix[i] = LeftEigenMatrix[i]/scale;

  // Lambda  = qn - can
  scale = 2.0*( Bn*(Bn*Bn - BdotB)/D - can*(UdotB - qn*Bn) );
  if ( abs(scale) < epsilon) scale = 1.0;
  theta1 = -Bn*Bcross/D - can*Ucross;
  LeftEigenMatrix[18] =  -( (-Ucross -can*Bcross/Bn)*(-w*can -Bz*Bn/D) + (w +can*Bz/Bn)*theta1) ;  
  LeftEigenMatrix[19] = -ny*( -w*can -Bz*Bn/D); LeftEigenMatrix[20] = nx*( -w*can -Bz*Bn/D);
  LeftEigenMatrix[21] = theta1; LeftEigenMatrix[22] = 0.0;
  LeftEigenMatrix[23] = -ny*(can*D/Bn)*(-w*can - Bz*Bn/D); LeftEigenMatrix[24] = nx*(can*D/Bn)*(-w*can - Bz*Bn/D);
  LeftEigenMatrix[25] = theta1*(can*D/Bn); LeftEigenMatrix[26] = 0.0;
  for(int i=18; i<=26; ++i) LeftEigenMatrix[i] = LeftEigenMatrix[i]/scale;

  // Lambda = qn - csn
  LeftEigenMatrix[28] = -u*(pow(csn,2.0) - pow(Bn,2.0)/D) + csn*(csn*csn*nx - Bn*Bx/D)/(gamma-1.0); LeftEigenMatrix[29] = -v*(pow(csn,2.0) - pow(Bn,2.0)/D) + csn*(csn*csn*ny - Bn*By/D)/(gamma-1.0);
  LeftEigenMatrix[30] = -w*(pow(csn,2.0) - pow(Bn,2.0)/D) + csn*(- Bn*Bz/D)/(gamma-1.0);
  LeftEigenMatrix[31] = (pow(csn,2.0) - pow(Bn,2.0)/D);
  LeftEigenMatrix[32] = -Bcross*ny*( ((gamma-2)/(gamma-1.0))*csn*csn - Bn*Bn/D); LeftEigenMatrix[33] = Bcross*nx*( ((gamma-2)/(gamma-1.0))*csn*csn - Bn*Bn/D);
  LeftEigenMatrix[34] = -Bz*( ((gamma-2)/(gamma-1.0))*csn*csn - Bn*Bn/D); LeftEigenMatrix[35] = 0.0;
  LeftEigenMatrix[27] = ( csn*csn/(gamma-1.0) + Bn*Bn/D -(DE+Pt)/D )*LeftEigenMatrix[31] -u*LeftEigenMatrix[28] -v*LeftEigenMatrix[29] -w*LeftEigenMatrix[30] -(Bx*LeftEigenMatrix[32] + By*LeftEigenMatrix[33] + Bz*LeftEigenMatrix[34])/D;
  scale = 0.0;
  for(int i=0; i<9; ++i) {
      scale += LeftEigenMatrix[27 + i] * RightEigenMatrix[3 + i*9];
  }
  if ( abs(scale) < epsilon ) scale = 1.0;
  for(int i=0;i<9;++i) LeftEigenMatrix[27 + i] = LeftEigenMatrix[27 + i]/scale;


  // Lambda = qn + csn
  LeftEigenMatrix[46] = -u*(pow(csn,2.0) - pow(Bn,2.0)/D) - csn*(csn*csn*nx - Bn*Bx/D)/(gamma-1.0); LeftEigenMatrix[47] = -v*(pow(csn,2.0) - pow(Bn,2.0)/D) - csn*(csn*csn*ny - Bn*By/D)/(gamma-1.0);
  LeftEigenMatrix[48] = -w*(pow(csn,2.0) - pow(Bn,2.0)/D) - csn*(- Bn*Bz/D)/(gamma-1.0);
  LeftEigenMatrix[49] = (pow(csn,2.0) - pow(Bn,2.0)/D);
  LeftEigenMatrix[50] = -Bcross*ny*( ((gamma-2)/(gamma-1.0))*csn*csn - Bn*Bn/D); LeftEigenMatrix[51] = Bcross*nx*( ((gamma-2)/(gamma-1.0))*csn*csn - Bn*Bn/D);
  LeftEigenMatrix[52] = -Bz*( ((gamma-2)/(gamma-1.0))*csn*csn - Bn*Bn/D); LeftEigenMatrix[53] = 0.0;
  LeftEigenMatrix[45] = ( csn*csn/(gamma-1.0) + Bn*Bn/D -(DE+Pt)/D )*LeftEigenMatrix[49] -u*LeftEigenMatrix[46] -v*LeftEigenMatrix[47] -w*LeftEigenMatrix[48] -(Bx*LeftEigenMatrix[50] + By*LeftEigenMatrix[51] + Bz*LeftEigenMatrix[52])/D;
  scale = 0.0;
  for(int i=0; i<9; ++i) {
      scale += LeftEigenMatrix[45 + i] * RightEigenMatrix[5 + i*9];
  }
  if ( abs(scale) < epsilon ) scale = 1.0;
  for(int i=0;i<9;++i) LeftEigenMatrix[45 + i] = LeftEigenMatrix[45 + i]/scale;

  // Lambda = qn - cfn
  LeftEigenMatrix[10] = -u*(pow(cfn,2.0) - pow(Bn,2.0)/D) + cfn*(cfn*cfn*nx - Bn*Bx/D)/(gamma-1.0); LeftEigenMatrix[11] = -v*(pow(cfn,2.0) - pow(Bn,2.0)/D) + cfn*(cfn*cfn*ny - Bn*By/D)/(gamma-1.0);
  LeftEigenMatrix[12] = -w*(pow(cfn,2.0) - pow(Bn,2.0)/D) + cfn*(- Bn*Bz/D)/(gamma-1.0);
  LeftEigenMatrix[13] = (pow(cfn,2.0) - pow(Bn,2.0)/D);
  LeftEigenMatrix[14] = -Bcross*ny*( ((gamma-2)/(gamma-1.0))*cfn*cfn - Bn*Bn/D); LeftEigenMatrix[15] = Bcross*nx*( ((gamma-2)/(gamma-1.0))*cfn*cfn - Bn*Bn/D);
  LeftEigenMatrix[16] = -Bz*( ((gamma-2)/(gamma-1.0))*cfn*cfn - Bn*Bn/D); LeftEigenMatrix[17] = 0.0;
  LeftEigenMatrix[9] = ( cfn*cfn/(gamma-1.0) + Bn*Bn/D -(DE+Pt)/D )*LeftEigenMatrix[13] -u*LeftEigenMatrix[10] -v*LeftEigenMatrix[11] -w*LeftEigenMatrix[12] -(Bx*LeftEigenMatrix[14] + By*LeftEigenMatrix[15] + Bz*LeftEigenMatrix[16])/D;
  scale = 0.0;
  for(int i=0; i<9; ++i) {
      scale += LeftEigenMatrix[9 + i] * RightEigenMatrix[1 + i*9];
  }
  if ( abs(scale) < epsilon ) scale = 1.0;
  for(int i=0;i<9;++i) LeftEigenMatrix[9 + i] = LeftEigenMatrix[9 + i]/scale;


  // Lambda = qn + cfn
  LeftEigenMatrix[64] = -u*(pow(cfn,2.0) - pow(Bn,2.0)/D) - cfn*(cfn*cfn*nx - Bn*Bx/D)/(gamma-1.0); LeftEigenMatrix[65] = -v*(pow(cfn,2.0) - pow(Bn,2.0)/D) - cfn*(cfn*cfn*ny - Bn*By/D)/(gamma-1.0);
  LeftEigenMatrix[66] = -w*(pow(cfn,2.0) - pow(Bn,2.0)/D) - cfn*(- Bn*Bz/D)/(gamma-1.0);
  LeftEigenMatrix[67] = (pow(cfn,2.0) - pow(Bn,2.0)/D);
  LeftEigenMatrix[68] = -Bcross*ny*( ((gamma-2)/(gamma-1.0))*cfn*cfn - Bn*Bn/D); LeftEigenMatrix[69] = Bcross*nx*( ((gamma-2)/(gamma-1.0))*cfn*cfn - Bn*Bn/D);
  LeftEigenMatrix[70] = -Bz*( ((gamma-2)/(gamma-1.0))*cfn*cfn - Bn*Bn/D); LeftEigenMatrix[71] = 0.0;
  LeftEigenMatrix[63] = ( cfn*cfn/(gamma-1.0) + Bn*Bn/D -(DE+Pt)/D )*LeftEigenMatrix[67] -u*LeftEigenMatrix[64] -v*LeftEigenMatrix[65] -w*LeftEigenMatrix[66] -(Bx*LeftEigenMatrix[68] + By*LeftEigenMatrix[69] + Bz*LeftEigenMatrix[70])/D;
  scale = 0.0;
  for(int i=0; i<9; ++i) {
      scale += LeftEigenMatrix[63 + i] * RightEigenMatrix[7 + i*9];
  }
  if ( abs(scale) < epsilon ) scale = 1.0;
  for(int i=0;i<9;++i) LeftEigenMatrix[63 + i] = LeftEigenMatrix[63 + i]/scale;

  // Lambda = ch
  LeftEigenMatrix[73] = pow(ch-qn,2.0)*((ch-qn)*nx/(gamma-1.0) - u) - (Bn/D)*((ch-gamma*qn)*nx/(gamma-1.0) - (By*u-Bx*v)*ny) ;
  LeftEigenMatrix[74] = pow(ch-qn,2.0)*((ch-qn)*ny/(gamma-1.0) - v) - (Bn/D)*((ch-gamma*qn)*ny/(gamma-1.0) - (By*u-Bx*v)*nx);
  LeftEigenMatrix[75] = -w*( pow(ch-qn,2.0) - Bn*Bn/D ) -Bz*Bn*(ch-qn)/(D*(gamma-1.0)) ;
  LeftEigenMatrix[76] = pow(ch - qn,2.0) - Bn*Bn/D;
  LeftEigenMatrix[77] = ch*nx*( pow(ch-qn,2.0) - Bn*Bn/D ) + Bcross*ny*(((2.0-gamma)/(gamma-1.0))*pow(ch-qn,2.0) + Bn*Bn/D );
  LeftEigenMatrix[78] = ch*ny*( pow(ch-qn,2.0) - Bn*Bn/D ) - Bcross*nx*(((2.0-gamma)/(gamma-1.0))*pow(ch-qn,2.0) + Bn*Bn/D );
  LeftEigenMatrix[79] = Bz *( ((2.0-gamma)/(gamma-1.0))*pow(ch-qn,2.0) + Bn*Bn/D ) ;
  LeftEigenMatrix[80] = pow(ch - qn,2.0) - Bn*Bn/D;
  LeftEigenMatrix[72] = (pow(ch-qn,2.0)/(gamma-1.0) + Bn*Bn/D -(DE+Pt)/D)*LeftEigenMatrix[76] + ch*(Bn/D)*LeftEigenMatrix[80] -u*LeftEigenMatrix[73] -v*LeftEigenMatrix[74] -w*LeftEigenMatrix[75] -(Bx*LeftEigenMatrix[77] + By*LeftEigenMatrix[78] + Bz*LeftEigenMatrix[79])/D  ;
  scale = 0.0;
  for(int i=0; i<9; ++i) {
      scale += LeftEigenMatrix[72 + i] * RightEigenMatrix[8 + i*9];
  }
  if ( abs(scale) < epsilon ) scale = 1.0;
  for(int i=0;i<9;++i) LeftEigenMatrix[72 + i] = LeftEigenMatrix[72 + i]/scale;

  // Lambda = -ch
  LeftEigenMatrix[1] = pow(-ch - qn,2.0)*((-ch - qn)*nx/(gamma-1.0) - u) - (Bn/D)*((-ch -gamma*qn)*nx/(gamma-1.0) - (By*u-Bx*v)*ny);
  LeftEigenMatrix[2] = pow(-ch - qn,2.0)*((-ch - qn)*ny/(gamma-1.0) - v) - (Bn/D)*((-ch -gamma*qn)*ny/(gamma-1.0) - (By*u-Bx*v)*nx);
  LeftEigenMatrix[3] = -w*( pow(-ch - qn,2.0) - Bn*Bn/D ) -Bz*Bn*(-ch - qn)/(D*(gamma-1.0)) ;
  LeftEigenMatrix[4] = pow(-ch - qn,2.0) - Bn*Bn/D;
  LeftEigenMatrix[5] = -ch*nx*( pow(-ch - qn,2.0) - Bn*Bn/D ) + Bcross*ny*(((2.0-gamma)/(gamma-1.0))*pow(-ch - qn,2.0) + Bn*Bn/D );
  LeftEigenMatrix[6] = -ch*ny*( pow(-ch - qn,2.0) - Bn*Bn/D ) - Bcross*nx*(((2.0-gamma)/(gamma-1.0))*pow(-ch - qn,2.0) + Bn*Bn/D );
  LeftEigenMatrix[7] = Bz *( ((2.0-gamma)/(gamma-1.0))*pow(-ch - qn,2.0) + Bn*Bn/D ) ;
  LeftEigenMatrix[8] = pow(-ch - qn,2.0) - Bn*Bn/D;
  LeftEigenMatrix[0] = (pow(-ch - qn,2.0)/(gamma-1.0) + Bn*Bn/D -(DE+Pt)/D)*LeftEigenMatrix[4] + -ch*(Bn/D)*LeftEigenMatrix[8] -u*LeftEigenMatrix[1] -v*LeftEigenMatrix[2] -w*LeftEigenMatrix[3] -(Bx*LeftEigenMatrix[5] + By*LeftEigenMatrix[6] + Bz*LeftEigenMatrix[7])/D  ;
  scale = 0.0;
  for(int i=0; i<9; ++i) {
      scale += LeftEigenMatrix[0 + i] * RightEigenMatrix[0 + i*9];
  }
  if ( abs(scale) < epsilon ) scale = 1.0;
  for(int i=0;i<9;++i) LeftEigenMatrix[0 + i] = LeftEigenMatrix[0 + i]/scale;

  
  return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to compute Characteristic Variables.
 * 
 * @Param array V Set of conservative variables..
 * @Param array C The Characteristic Variable.
 * @Param I Identifier for Characteristic Variable.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::convertVariabletoCharacteristic(int *V, int *C, unsigned I) {
    double Sum = 0.0;
    for(int k=0; k < Dimension; ++k) {
        double *VarK = variable[C[k]];
        for(int i=0; i< (N+1)*(N+1); ++i) {
        Sum = 0.0;
        for(int j=0; j < Dimension; ++j){
            Sum += variable[V[j]][i]*LeftEigenMatrix[k*Dimension +j];
        }     
        VarK[i] = Sum;
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
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::convertCharacteristictoVariable(int *C, int *V, unsigned I) {
    double Sum = 0.0;
    for(int k=0; k < Dimension; ++k) {
        double *VarK = variable[V[k]];
        for(int i=0; i< (N+1)*(N+1); ++i) {
         Sum = 0.0;
         for(int j=0; j < Dimension; ++j){
            Sum += variable[C[j]][i]*RightEigenMatrix[k*Dimension +j];
          }     
         VarK[i] = Sum;
        } 
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

        double *BRv = boundaryRight[v], *NRv = neighboringRight[v];
        double *BLv = boundaryLeft[v], *NLv = neighboringLeft[v];                                                     
        for(int i=0; i<=N; i++){
            numericalFlux[i*(N+1)+N] = 0.5*(BRv[i*(N+1)] + NRv[i*(N+1)] );   
            numericalFlux[i*(N+1)]   = 0.5*(BLv[i*(N+1)]  + NLv[i*(N+1)]);   
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

        double *BRv = boundaryRight[v], *NRv = neighboringRight[v], *BRfv = boundaryRight[fluxVariable], *NRfv = neighboringRight[fluxVariable], *BRcv = boundaryRight[conserVar], *NRcv = neighboringRight[conserVar];
        double *BLv = boundaryLeft[v], *NLv = neighboringLeft[v], *BLfv = boundaryLeft[fluxVariable], *NLfv = neighboringLeft[fluxVariable], *BLcv = boundaryLeft[conserVar], *NLcv = neighboringLeft[conserVar];                                                     
        for(int i=0; i<=N; i++){
            numericalFlux[i*(N+1)+N] = 0.5*(BRv[i*(N+1)] + NRv[i*(N+1)] + MAX(fabs(BRfv[i*(N+1)]), fabs(NRfv[i*(N+1)]))*(BRcv[i*(N+1)] - NRcv[i*(N+1)])  ) ;   
            numericalFlux[i*(N+1)]   = 0.5*(BLv[i*(N+1)]  + NLv[i*(N+1)]  - MAX(fabs(BLfv[i*(N+1)]), fabs(NLfv[i*(N+1)]))*(BLcv[i*(N+1)] - NLcv[i*(N+1)])  ) ;   
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

        double *BTv = boundaryTop[v], *NTv = neighboringTop[v] ;
        double *BBv = boundaryBottom[v], *NBv = neighboringBottom[v];                                                                                                          
        for(int i=0; i<=N; i++){
            numericalFlux[N*(N+1)+i]= 0.5*(BTv[i] + NTv[i]);
            numericalFlux[i]        = 0.5*(BBv[i] + NBv[i]); 
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

        double *BTv = boundaryTop[v], *NTv = neighboringTop[v], *BTfv = boundaryTop[fluxVariable], *NTfv = neighboringTop[fluxVariable], *BTcv = boundaryTop[conserVar], *NTcv = neighboringTop[conserVar];
        double *BBv = boundaryBottom[v], *NBv = neighboringBottom[v], *BBfv = boundaryBottom[fluxVariable], *NBfv = neighboringBottom[fluxVariable], *BBcv = boundaryBottom[conserVar], *NBcv = neighboringBottom[conserVar];                                                                                                          
        for(int i=0; i<=N; i++){
            numericalFlux[N*(N+1)+i]= 0.5*(BTv[i] + NTv[i] + MAX(fabs(BTfv[i]), fabs(NTfv[i]))*(BTcv[i] - NTcv[i]));
            numericalFlux[i]        = 0.5*(BBv[i] + NBv[i] - MAX(fabs(BBfv[i]), fabs(NBfv[i]))*(BBcv[i] - NBcv[i])); 
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
 * @Param conserVar Corresponding Conservative variable.
 * @Param fluxType  The type of flux that is to be used. eg "central"
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::delByDelX(int *V, int *VDash, int *ConserVar, string fluxType, int *FluxVariable, unsigned size) {
    double dy = (y_end - y_start);
    double dx = (x_end - x_start);
    
    if(fluxType == "central") {
        int v, vDash, conserVar, fluxVariable = FluxVariable[0];
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        
        for(int temp =0 ; temp < size; ++temp) {
            v = V[temp]; vDash = VDash[temp]; conserVar = ConserVar[temp];
            zeros(numericalFlux, (N+1)*(N+1)); 

            double *BRv = boundaryRight[v], *NRv = neighboringRight[v];
            double *BLv = boundaryLeft[v], *NLv = neighboringLeft[v];                                                     
            for(int i=0; i<=N; i++){
               numericalFlux[i*(N+1)+N] = 0.5*(BRv[i*(N+1)] + NRv[i*(N+1)]);   
               numericalFlux[i*(N+1)]   = 0.5*(BLv[i*(N+1)]  + NLv[i*(N+1)]);   
            }
            /// vDash = -0.5*dy*D*v
            cblas_dgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dy, derivativeMatrix_x, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

            /// Adding the numeical Flux terms as necessary.
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dy, fluxMatrix_right,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dy, fluxMatrix_left,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

            /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);
        }
        
        delete[] numericalFlux;
        delete[] auxillaryVariable;
    }

    else if(fluxType == "rusanov") {
        int v, vDash, conserVar, fluxVariable = FluxVariable[0];
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        
        for(int temp =0 ; temp < size; ++temp) {
            v = V[temp]; vDash = VDash[temp]; conserVar = ConserVar[temp];
            zeros(numericalFlux, (N+1)*(N+1)); 

            double *BRv = boundaryRight[v], *NRv = neighboringRight[v], *BRfv = boundaryRight[fluxVariable], *NRfv = neighboringRight[fluxVariable], *BRcv = boundaryRight[conserVar], *NRcv = neighboringRight[conserVar];
            double *BLv = boundaryLeft[v], *NLv = neighboringLeft[v], *BLfv = boundaryLeft[fluxVariable], *NLfv = neighboringLeft[fluxVariable], *BLcv = boundaryLeft[conserVar], *NLcv = neighboringLeft[conserVar];                                                     
            for(int i=0; i<=N; i++){
               numericalFlux[i*(N+1)+N] = 0.5*(BRv[i*(N+1)] + NRv[i*(N+1)] + MAX(fabs(BRfv[i*(N+1)]), fabs(NRfv[i*(N+1)]))*(BRcv[i*(N+1)] - NRcv[i*(N+1)])  ) ;   
               numericalFlux[i*(N+1)]   = 0.5*(BLv[i*(N+1)]  + NLv[i*(N+1)]  - MAX(fabs(BLfv[i*(N+1)]), fabs(NLfv[i*(N+1)]))*(BLcv[i*(N+1)] - NLcv[i*(N+1)])  ) ;   
            }

            /// vDash = -0.5*dy*D*v
            cblas_dgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dy, derivativeMatrix_x, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

            /// Adding the numeical Flux terms as necessary.
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dy, fluxMatrix_right,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dy, fluxMatrix_left,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

            /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);
        }
    
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
void DG_Element_2d::delByDelY(int *V, int *VDash, int *ConserVar, string fluxType, int *FluxVariable, unsigned size) {
    double dy = (y_end - y_start);
    double dx = (x_end - x_start);

    if(fluxType == "central") {
        int v, vDash, conserVar, fluxVariable = FluxVariable[0];
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        
        for(int temp =0 ; temp < size; ++temp) {
            v = V[temp]; vDash = VDash[temp]; conserVar = ConserVar[temp];
            zeros(numericalFlux, (N+1)*(N+1));                                                       
            double *BTv = boundaryTop[v], *NTv = neighboringTop[v];
            double *BBv = boundaryBottom[v], *NBv = neighboringBottom[v];                                                                                                          
            for(int i=0; i<=N; i++){
               numericalFlux[N*(N+1)+i]= 0.5*(BTv[i] + NTv[i]);
               numericalFlux[i]        = 0.5*(BBv[i] + NBv[i]); 
            }
            /// vDash = -0.5*dy*D*v
            cblas_dgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dx, derivativeMatrix_y, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

            /// Adding the numeical Flux terms as necessary.
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dx, fluxMatrix_top,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dx, fluxMatrix_bottom,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

            /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);
        }
        
        delete[] numericalFlux;
        delete[] auxillaryVariable;
    }
    
    else if(fluxType == "rusanov") {
        int v, vDash, conserVar, fluxVariable = FluxVariable[0];
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        
        for(int temp =0 ; temp < size; ++temp) {
            v = V[temp]; vDash = VDash[temp]; conserVar = ConserVar[temp];
            zeros(numericalFlux, (N+1)*(N+1));                                                       
            double *BTv = boundaryTop[v], *NTv = neighboringTop[v], *BTfv = boundaryTop[fluxVariable], *NTfv = neighboringTop[fluxVariable], *BTcv = boundaryTop[conserVar], *NTcv = neighboringTop[conserVar];
            double *BBv = boundaryBottom[v], *NBv = neighboringBottom[v], *BBfv = boundaryBottom[fluxVariable], *NBfv = neighboringBottom[fluxVariable], *BBcv = boundaryBottom[conserVar], *NBcv = neighboringBottom[conserVar];                                                                                                          
            for(int i=0; i<=N; i++){
               numericalFlux[N*(N+1)+i]= 0.5*(BTv[i] + NTv[i] + MAX(fabs(BTfv[i]), fabs(NTfv[i]))*(BTcv[i] - NTcv[i]));
               numericalFlux[i]        = 0.5*(BBv[i] + NBv[i] - MAX(fabs(BBfv[i]), fabs(NBfv[i]))*(BBcv[i] - NBcv[i])); 
            }
            /// vDash = -0.5*dy*D*v
            cblas_dgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dx, derivativeMatrix_y, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

            /// Adding the numeical Flux terms as necessary.
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dx, fluxMatrix_top,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dx, fluxMatrix_bottom,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

            /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
            cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);

        }
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

        double *BRv = boundaryRight[v], *NRv = neighboringRight[v];
        double *BLv = boundaryLeft[v], *NLv = neighboringLeft[v];                                                     
        for(int i=0; i<=N; i++){
            numericalFlux[i*(N+1)+N] = 0.5*(BRv[i*(N+1)] + NRv[i*(N+1)] );   
            numericalFlux[i*(N+1)]   = 0.5*(BLv[i*(N+1)]  + NLv[i*(N+1)]);   
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
        
        double *BTv = boundaryTop[v], *NTv = neighboringTop[v] ;
        double *BBv = boundaryBottom[v], *NBv = neighboringBottom[v];                                                                                                          
        for(int i=0; i<=N; i++){
            numericalFlux[N*(N+1)+i]= 0.5*(BTv[i] + NTv[i]);
            numericalFlux[i]        = 0.5*(BBv[i] + NBv[i]); 
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