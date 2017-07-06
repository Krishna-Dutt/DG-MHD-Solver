#include "../../includes/DG_Element_2d/DG_Element_2d.h"
#include "../../includes/Utilities/Zeros.h"

#include "../../includes/Utilities/Inverse.h"

#include "../../includes/Utilities/FluxMatrix.h"
#include "../../includes/Utilities/MassMatrix.h"
#include "../../includes/Utilities/DerivativeMatrix.h"
#include "../../includes/Utilities/LobattoNodes.h"
#include "../../includes/Utilities/MinMod.h"

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

    for ( map<string, double*>::iterator itr = variable.begin() ;itr != variable.end(); itr++){
      delete[] (itr->second);
    }
    // Do I need to delete other maps pointing to Boundary elements, since no new memory is dynamically allocated for them ??
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
 * @Synopsis  This functions creates space in order to take in one more variable on which operators are needed to be
 * applied.
 *
 * @Param v  This is the name of the variable which is to be added.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::addVariable_withoutBoundary(string v) {
    double * newVariable = new double[(N+1)*(N+1)]; /// Allocating the space for the new variable which is to be created.
    variable[v] = newVariable; /// Now assigning the same to the map.

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
void DG_Element_2d::addVariable_CellCentered(string v) {
    double * newVariable = new double ;/// Allocating the space for the new variable which is to be created.
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
void DG_Element_2d::ResetVariables_CellCentered(string v, double value) {
   *variable[v] = value; 
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function resets the Map to the Outflow Boundaries.
 *
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::ResetMap_OutFlow() {
    OutFlow["Top"] = false ;
    OutFlow["Bottom"] = false;
    OutFlow["Right"]= false;
    OutFlow["Left"] = false;
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
void DG_Element_2d::updateOutFlowBoundary(string u, string v) {
  double start = -1.0;
  double end = 1.0;
  
  if (lobattoIntegration(start, end, N, boundaryTop[v]) < 0.0) {
    OutFlow["Top"] = true;
  }
  if (lobattoIntegration(start, end, N, boundaryBottom[v]) > 0.0) {
    OutFlow["Bottom"] = true;
  }
  
  if (lobattoIntegration(start, end, N, boundaryLeft[u]) > 0.0) {
    OutFlow["Left"] = true;
  }
  if (lobattoIntegration(start, end, N, boundaryRight[u]) < 0.0) {
    OutFlow["Right"] = true;
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
void DG_Element_2d::updateCellMarker(string v, string m) {
  double radius = 1.0;
  double OutflowSize = 0.0;
  double MaxVariable = 0.0;
  double VariableFlux = 0.0;

  if (OutFlow["Top"]) {
    VariableFlux += lobattoIntegration(x_start, x_end, N, boundaryTop[v]);
    VariableFlux -= lobattoIntegration(x_start, x_end, N, neighboringTop[v]);
    OutflowSize += abs(x_end -x_start);

    for(int i=0; i<=N; ++i) {
        MaxVariable = MAX(MaxVariable, *boundaryTop[v][i]);
    }
  }
  if (OutFlow["Bottom"]) {
    VariableFlux += lobattoIntegration(x_start, x_end, N, boundaryBottom[v]);
    VariableFlux -= lobattoIntegration(x_start, x_end, N, neighboringBottom[v]);
    OutflowSize += abs(x_end -x_start);

    for(int i=0; i<=N; ++i) {
        MaxVariable = MAX(MaxVariable, *boundaryBottom[v][i]);
    }
  }
  if (OutFlow["Left"]) {
    VariableFlux += lobattoIntegration(y_start, y_end, N, boundaryLeft[v]);
    VariableFlux -= lobattoIntegration(y_start, y_end, N, neighboringLeft[v]);
    OutflowSize += abs(y_end -y_start);

    for(int i=0; i<=N; ++i) {
        MaxVariable = MAX(MaxVariable, *boundaryLeft[v][i]);
    }
  }
  if (OutFlow["Right"]) {
    VariableFlux += lobattoIntegration(y_start, y_end, N, boundaryRight[v]);
    VariableFlux -= lobattoIntegration(y_start, y_end, N, neighboringRight[v]);
    OutflowSize += abs(y_end -y_start);

    for(int i=0; i<=N; ++i) {
        MaxVariable = MAX(MaxVariable, *boundaryRight[v][i]);
    }
  }

  MaxVariable = variable[v][0];
  for (int i = 0; i < (N+1)*(N+1) ; ++i) {
    MaxVariable = MAX(MaxVariable,variable[v][i]);
  }


// Assert that MaxVariable is never equal to ZERO !!  

  radius = MIN(abs(x_start-x_end),abs(y_start-y_end)) * 0.5;

  *variable[m] = abs(VariableFlux) / ( abs(OutflowSize) * MaxVariable * pow(radius, 0.5 * (N+1)));
  if (*variable[m] > 1.0) {
    *variable[m] = 1.0;
  }
  else {
    *variable[m] = 0.0;
  }
  
  for(int b=0; b <(N+1)*(N+1); ++b) {
      variable["CellMarkerGlobal"][b] = *variable[m];
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
void DG_Element_2d::computeMoments(string v, string m) {
  /// Multiplying inverse of VanderMand Matrix with the variable array to get the corresponding moments.
  cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 1.0, inverseVanderMandMatrix,(N+1)*(N+1), variable[v],1,0,variable[m],1);

// Additional filter to set values less than 1e-8 to zero !!, Recheck /Temporary
/*
  for (int i=0; i <(N+1)*(N+1); ++i) {
      if (abs(variable[m][i]) <= 1e-7) {
          variable[m][i] = 0.0;
      }
  }
*/
  return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function computes the nodals values of the variable given its moments
 * using VandeMandMatrix with Legendre Basis.
 *
 * @Param m This gives the moments of the variable
 * @Param v This is the variable whose nodal values are to be computed.
 * @Param cm This is the cell marker used to identified troubled cells.
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::convertMomentToVariable(string m, string v, string cm) {
  /// Multiplying  VanderMand Matrix with the moments to obtained the nodal values of the variable.

 //if (*variable[cm])
  { // Checking if cell marker is not equal to zero
  //cout << "Calling :: convertMomentToVariable()\n";
  
  
  cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 1.0, vanderMandMatrix,(N+1)*(N+1), variable[m],1,0,variable[v],1);
  

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
*/
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::limitMoments(string m, string modm, string cm) {

//if (*variable[cm]) 
  { // Checking if cell marker is not equal to zero
    int count, Tempi, Tempj, i, j;
    count = N+1;
    double Temp1, Temp2, AlphaN;
    //AlphaN = sqrt((2.0*N)/(2.0*N+1));
    AlphaN = sqrt((2.0*N -1.0)/(2.0*N +1));  // Similar to a diffusion coefficient
    //AlphaN = 0.5/sqrt(4.0*N*N - 1.0);        // Too diffusive
    //AlphaN = (2.0*N-0.5)/sqrt(4.0*N*N-1.0);  // Too diffusive
    //AlphaN = 1.0;                           // Do not work for HO polynomials
    //AlphaN = 1.0/(2.0*N-1.0);               // Do not work for HO polynomials

    for(i=(N+1)*(N+1)-1; i > 0; i = i - (N+2)) {
     --count;
     AlphaN = sqrt((2.0*(count)-1.0)/(2.0*(count)+1.0));
     for(j=0; j < count; ++j) {
       Tempi = i-j;
       Tempj = i - j*(N+1);
       // Original minmod detector
       Temp1 = MinMod(variable[m][Tempi], AlphaN*(rightNeighbor->variable[m][Tempi-1] -variable[m][Tempi-1]), AlphaN*(variable[m][Tempi-1] -leftNeighbor->variable[m][Tempi-1]) , AlphaN*(topNeighbor->variable[m][Tempi-(N+1)] -variable[m][Tempi-(N+1)]), AlphaN*(variable[m][Tempi-(N+1)] -bottomNeighbor->variable[m][Tempi-(N+1)]));
       Temp2 = MinMod(variable[m][Tempj], AlphaN*(rightNeighbor->variable[m][Tempj-1] -variable[m][Tempj-1]), AlphaN*(variable[m][Tempj-1] -leftNeighbor->variable[m][Tempj-1]) , AlphaN*(topNeighbor->variable[m][Tempj-(N+1)] -variable[m][Tempj-(N+1)]), AlphaN*(variable[m][Tempj-(N+1)] -bottomNeighbor->variable[m][Tempj-(N+1)]));
       
       if ( Temp1 != variable[m][Tempi] || Temp2 != variable[m][Tempj] ) {
         variable[modm][Tempi] = Temp1;
         variable[modm][Tempj] = Temp2;
       }
       else //if( Temp1 !=0 && Temp2 !=0) 
       {
         return ; // Need to exit both loops
       }
     }
     // Special Case for end values, when Tempi or Tempj access zero order polynomials !!
       Tempi = i-j;
       Tempj = i - j*(N+1);
       
       // Original Detector
       Temp1 = MinMod(variable[m][Tempi], AlphaN*(topNeighbor->variable[m][Tempi-(N+1)] -variable[m][Tempi-(N+1)]), AlphaN*(variable[m][Tempi-(N+1)] -bottomNeighbor->variable[m][Tempi-(N+1)]));
       Temp2 = MinMod(variable[m][Tempj], AlphaN*(rightNeighbor->variable[m][Tempj-1] -variable[m][Tempj-1]), AlphaN*(variable[m][Tempj-1] -leftNeighbor->variable[m][Tempj-1]));
       
       if ( Temp1 != variable[m][Tempi] || Temp2 != variable[m][Tempj] ) {
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
 * @Synopsis  This functions creates space in order to take in one more variable on which operators are needed to be
 * applied.
 *
 * @Param v  This is the name of the variable which is to be added.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::addVariable_withBoundary(string v) {
    double * newVariable = new double[(N+1)*(N+1)]; /// Allocating the space for the new variable which is to be created.
    variable[v] = newVariable; /// Now assigning the same to the map.
    
    // **b_top is used because it will store the address to the boundary element, So whenever the actual value in the double* of the variable is changed then this will also change automatically. The same holds for all the other following mentioned variables.
    
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
        n_top[i]    =   b_top[i]    = &(variable[v][N*(N+1)+i]);
        n_left[i]   =   b_left[i]   = &(variable[v][i*(N+1)]);
        n_right[i]  =   b_right[i]  = &(variable[v][i*(N+1)+N]);
    }


    boundaryTop[v]      = b_top;
    boundaryRight[v]    = b_right;
    boundaryBottom[v]   = b_bottom;
    boundaryLeft[v]     = b_left;


    neighboringTop[v]   = n_top;
    neighboringRight[v] = n_right;
    neighboringBottom[v]= n_bottom;
    neighboringLeft[v]  = n_left;

    boundaryVariables.push_back(v);

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
void DG_Element_2d::addVariable_onlyBoundary(string v) {
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
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to set the value of a variable with the help of a function.
 *
 * @Param v This is the name of the variable whose value is to be  set.
 * @Param f This is the f(x, y) which is used to initialize the variable.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::initializeVariable(string v, function<double(double, double)> f) {
    for(int i=0; i<((N+1)*(N+1)); i++)
        variable[v][i] = f(X[i], Y[i]);
    
    return ;
}


void DG_Element_2d::setVariableNeighbors(string v) {
    int j;
    if(topNeighbor!=this) {
        for( j = 0 ; j <= N; j++) 
            neighboringTop[v][j] = topNeighbor->boundaryBottom[v][j];
    }
    if(rightNeighbor!=this) {
        for( j = 0 ; j <= N; j++) 
            neighboringRight[v][j] = rightNeighbor->boundaryLeft[v][j];
    }
    if(leftNeighbor!=this) {
        for( j = 0 ; j <= N; j++) 
            neighboringLeft[v][j] = leftNeighbor->boundaryRight[v][j];
    }
    if(bottomNeighbor!=this) {
        for( j = 0 ; j <= N; j++) 
            neighboringBottom[v][j] = bottomNeighbor->boundaryTop[v][j];
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
    string currentVariable;
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
 * @Param fluxType  The type of flux that is to be used. eg "central"
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::delByDelX(string v, string vDash, string conserVar, string fluxType, string fluxVariable = "") {
    double dy = (y_end - y_start);
    double dx = (x_end - x_start);
    
    if(fluxType == "central") {
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        zeros(numericalFlux, (N+1)*(N+1));                                                       
        for(int i=0; i<=N; i++){
            numericalFlux[i*(N+1)+N]    = 0.5*( *boundaryRight[v][i]    + *neighboringRight[v][i] ) ;   
            numericalFlux[i*(N+1)]    = 0.5*( *boundaryLeft[v][i]     + *neighboringLeft[v][i] ) ;  
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
            numericalFlux[i*(N+1)+N] = 0.5*(*boundaryRight[v][i] + *neighboringRight[v][i] + MAX(fabs(*boundaryRight[fluxVariable][i]), fabs(*neighboringRight[fluxVariable][i]))*(*boundaryRight[conserVar][i] - *neighboringRight[conserVar][i])  ) ;   
            numericalFlux[i*(N+1)]   = 0.5*(*boundaryLeft[v][i]  + *neighboringLeft[v][i]  - MAX(fabs(*boundaryLeft[fluxVariable][i]), fabs(*neighboringLeft[fluxVariable][i]))*(*boundaryLeft[conserVar][i] - *neighboringLeft[conserVar][i])  ) ;   
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
void DG_Element_2d::delByDelY(string v, string vDash, string conserVar, string fluxType, string fluxVariable = "") {
    double dy = (y_end - y_start);
    double dx = (x_end - x_start);

    if(fluxType == "central") {
        double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
        double* auxillaryVariable    =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        zeros(numericalFlux, (N+1)*(N+1));                                                       
        for(int i=0; i<=N; i++){
            numericalFlux[N*(N+1)+i]    = 0.5*( *boundaryTop[v][i]    + *neighboringTop[v][i] ) ;   
            numericalFlux[i]            = 0.5*( *boundaryBottom[v][i]     + *neighboringBottom[v][i] ) ;  
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
            numericalFlux[N*(N+1)+i]= 0.5*(*boundaryTop[v][i] + *neighboringTop[v][i] + MAX(fabs(*boundaryTop[fluxVariable][i]), fabs(*neighboringTop[fluxVariable][i]))*(*boundaryTop[conserVar][i] - *neighboringTop[conserVar][i]));
            numericalFlux[i]        = 0.5*(*boundaryBottom[v][i] + *neighboringBottom[v][i] - MAX(fabs(*boundaryBottom[fluxVariable][i]), fabs(*neighboringBottom[fluxVariable][i]))*(*boundaryBottom[conserVar][i] - *neighboringBottom[conserVar][i])); 
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
void DG_Element_2d::axpy(double a, string x, string y) {
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
void DG_Element_2d::scal(double a, string x) {
    cblas_dscal((N+1)*(N+1), a, variable[x], 1);
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
void DG_Element_2d::setFunctionsForVariables(string x, string y, function<double(double, double)> f, string z) {
    for(int i = 0 ; i < (N+1)*(N+1); i++)
        variable[z][i] = f(variable[x][i],variable[y][i]);
    return ;
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
void DG_Element_2d::setFunctionsForVariables(string w, string x, string y, function<double(double, double, double)> f, string z) {
    for(int i = 0 ; i < (N+1)*(N+1); i++)
        variable[z][i] = f(variable[w][i],variable[x][i],variable[y][i]);
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of Boundary variable z to f(x, y).
 * 
 *
 * @Param x The first parameter of the function.
 * @Param y The second parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setFunctionsForBoundaryVariables(string x, string y, function<double(double, double)> f, string z) {
    for(int i = 0 ; i <= N ; i++){
         *boundaryTop[z][i] = f(*boundaryTop[x][i], *boundaryTop[y][i]);
         *boundaryRight[z][i] = f(*boundaryRight[x][i], *boundaryRight[y][i]);
         *boundaryLeft[z][i] = f(*boundaryLeft[x][i], *boundaryLeft[y][i]);
         *boundaryBottom[z][i] = f(*boundaryBottom[x][i], *boundaryBottom[y][i]);
    }
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function used to change the value of  Boundary variable z to f(w, x, y).
 *
 * @Param w The first parameter of the function
 * @Param x The second parameter of the function.
 * @Param y The third parameter of the function.
 * @Param functionf The function `f` which is required for the intended mapping.
 * @Param z The variable in which the value is to be stored
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setFunctionsForBoundaryVariables(string w, string x, string y, function<double(double, double, double)> f, string z) {
    for(int i = 0 ; i <= N ; i++){
         *boundaryTop[z][i] = f(*boundaryTop[w][i], *boundaryTop[x][i], *boundaryTop[y][i]);
         *boundaryRight[z][i] = f(*boundaryRight[w][i], *boundaryRight[x][i], *boundaryRight[y][i]);
         *boundaryLeft[z][i] = f(*boundaryLeft[w][i], *boundaryLeft[x][i], *boundaryLeft[y][i]);
         *boundaryBottom[z][i] = f(*boundaryBottom[w][i], *boundaryBottom[x][i], *boundaryBottom[y][i]);
    }
    return;
}


double DG_Element_2d::l2Norm(string v1, string v2) {
    double* diff = new double[(N+1)*(N+1)];
    cblas_dscal((N+1)*(N+1), 0.0, diff, 1);
    cblas_daxpy((N+1)*(N+1),  1.0, variable[v1], 1, diff, 1);
    cblas_daxpy((N+1)*(N+1), -1.0, variable[v2], 1, diff, 1);
    double norm2 =  (cblas_dnrm2((N+1)*(N+1),diff, 1));    
    delete[] diff;
    return norm2;
}

/* ----------------------------------------------------------------------------*/
//--------------------------VIRTUAL FUNCTIONS----------------------------------//
/*-----------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates maps from a variable to its corresponding Boundary Setting methods
 *
 * @Param v  This is the name of the variable whose map is to be set.
 * @Param type This is the type of the boundary.
 * @Param b This defines the location of the boundary element.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::assignBoundary(string v, string type, char b) {
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates the linear system for a variable corresponding 
 * to its Boundary settings.
 *
 * @Param v  This is the name of the variable whose Boundary is to be set.
 * @Param b This defines the location of the boundary element.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setBoundaryValue(string v, string b) {
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function implements Dirichlet BC.
 *
 * @Param Matrix  RHS Matrix that is to be modified.
 * @Param List I List of indices and corresponding increments to access elements to be modified.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::DirichletBoundary(double *Matrix, initializer_list<int> I) {
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function implements Neumann BC ( To be precise, zero gradient condition).
 *
 * @Param Matrix  RHS Matrix that is to be modified.
 * @Param List I List of indices and corresponding increments required to access elements to be modified.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::NeumannBoundary(double *Matrix, initializer_list<int> I) {
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function implements Periodic BC .
 *
 * @Param Matrix  RHS Matrix that is to be modified.
 * @Param List I List of indices and corresponding increments required to access elements to be modified.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::PeriodicBoundary(double *Matrix, initializer_list<int> I) {
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates Neumann Boundaries of a cell.
 *
 * @Param v This is the variable ,whose boundary values are to be fixed.
 * @Param Matrix  RHS Matrix that is to be modified.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::updateNeumann(string v, double *Matrix) {
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates Dirichlet Boundaries of a cell.
 *
 * @Param v This is the variable ,whose boundary values are to be fixed.
 * @Param Matrix  RHS Matrix that is to be modified.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::updateDirichlet(string v, double *Matrix) {
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates Variables at the Boundary.
 *
 * @Param v This is the variable ,whose boundary values are to be fixed.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::updateBoundaryVariables(string v) {
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function performs MinMod for Boundary elements, eliminating all non-existing 
 * neighbouring cells.
 *
 * @Param m This string represents the moment of the variable
 * @Param Index This is represents to be limited.
 * @Param Alpha This is the scaling factor used in Lilia's Moment Limiter.
 * @Param R This is the pointer to RightNeighbour.
 * @Param L This is the pointer to LeftNeighbour.
 * @Param T This is the pointer to TopNeighbour.
 * @Param B This is the pointer to BottomNeighbour.
*/
/* ----------------------------------------------------------------------------*/
double DG_Element_2d::BoundaryMinMod(string m, int Index, double Alpha, DG_Element_2d* R, DG_Element_2d* L, DG_Element_2d* T, DG_Element_2d* B) {
    return 0.0 ;
}
