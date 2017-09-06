#include "../../includes/DG_BoundaryElement_2d/DG_BoundaryElement_2d.h"
#include "../../includes/Utilities/Zeros.h"
#include "../../includes/Utilities/Inverse.h"
#include "../../includes/Utilities/FluxMatrix.h"
#include "../../includes/Utilities/MassMatrix.h"
#include "../../includes/Utilities/DerivativeMatrix.h"
#include "../../includes/Utilities/LobattoNodes.h"
#include "../../includes/Utilities/MinMod.h"

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
DG_BoundaryElement_2d::DG_BoundaryElement_2d(int _N, double x1, double y1, double x2, double y2) : DG_Element_2d(_N, x1, y1, x2, y2) {

}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the destructor which deallocates the member variables of the class and destroys the object.
 *
*/
/* ----------------------------------------------------------------------------*/
DG_BoundaryElement_2d::~DG_BoundaryElement_2d() {
    for ( map<string, double*>::iterator itr = DirichletTop.begin() ;itr != DirichletTop.end(); itr++){
      delete[] (itr->second);
    }
    for ( map<string, double*>::iterator itr = DirichletBottom.begin() ;itr != DirichletBottom.end(); itr++){
      delete[] (itr->second);
    }
    for ( map<string, double*>::iterator itr = DirichletRight.begin() ;itr != DirichletRight.end(); itr++){
      delete[] (itr->second);
    }
    for ( map<string, double*>::iterator itr = DirichletLeft.begin() ;itr != DirichletLeft.end(); itr++){
      delete[] (itr->second);
    }

}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates maps from a variable to its corresponding Boundary Setting methods
 *
 * @Param v  This is the name of the variable whose map is to be set.
 * @Param type This is the type of the boundary.
 * @Param b This defines the location of the boundary element.
 */
/* ----------------------------------------------------------------------------*/
void DG_BoundaryElement_2d::assignBoundary(string v, string type, char b) {
     switch(b) {
        case 't' : // `t` or `T` for top 
        case 'T' :
            TopBoundary[v] = type;
            if (type == "dirichlet") {
                DirichletTop[v] = new double [N+1];
                //DirichletTop[v].resize(N+1);
                //DirichletTop[v] = type;
                for(int i=0; i<=N; i++){
                    DirichletTop[v][i] = *boundaryTop[v][i];
                }
            }
            break;
        case 'r' : // `r` or `R` for right
        case 'R' :
            RightBoundary[v] = type;
            if (type == "dirichlet") {
                DirichletRight[v] = new double [N+1];
                //DirichletRight[v].resize(N+1);
                 //DirichletRight[v] = type;
                for(int i=0; i<=N; i++){
                    DirichletRight[v][i] = *boundaryRight[v][i];
                }
            }
            break;
        case 'b' : // `b` or `B` for bottom
        case 'B' :
            BottomBoundary[v] = type;
            if (type == "dirichlet") {
                DirichletBottom[v] = new double [N+1];
                //DirichletBottom[v].resize(N+1);
                //DirichletBottom[v] = type;
                for(int i=0; i<=N; i++){
                    DirichletBottom[v][i] = *boundaryBottom[v][i];
                }
            }
            break;
        case 'l' : // `l` or `L` for left
        case 'L' :
            LeftBoundary[v] = type;
            if (type == "dirichlet") {
                DirichletLeft[v] = new double [N+1];
               //DirichletLeft[v].resize(N+1);
               //DirichletLeft[v] = type;
                for(int i=0; i<=N; i++){
                    DirichletLeft[v][i] = *boundaryLeft[v][i];
                }
            }
            break;
        default:
            cout << "WARNING!. No such neighbor type " << type << endl;
    }

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
void DG_BoundaryElement_2d::setBoundaryValue(string v, string b) {
    return ;
}
/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function implements Dirichlet BC.
 *
 * @Param Matrix  RHS Matrix that is to be modified.
 * @Param List I List of indices and corresponding increments to access elements to be modified, which 
 * includes the starting index, and the increment in index required.
 */
/* ----------------------------------------------------------------------------*/
void DG_BoundaryElement_2d::DirichletBoundary(double *Matrix, initializer_list<int> I) {
    int Start_I, Inc;
    Start_I = I.begin()[0];
    Inc = I.begin()[1];

    for (int i=0; i<=N; ++i) {
        for (int j=0; j<(N+1)*(N+1); ++j) {
            Matrix[(Start_I + i*Inc)*(N+1)*(N+1) + j] = 0.0;
        }
    }

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function implements Neumann BC ( To be precise, zero gradient condition).
 *
 * @Param Matrix  RHS Matrix that is to be modified.
 * @Param List I List of indices and corresponding increments required to access elements to be modified, which 
 * includes the starting index, the increment in index and the  index from which to copy data.
 */
/* ----------------------------------------------------------------------------*/
void DG_BoundaryElement_2d::NeumannBoundary(double *Matrix, initializer_list<int> I) {
    int Start_I, Inc, Copy_I;
    Start_I = I.begin()[0];
    Inc = I.begin()[1];
    Copy_I = I.begin()[2];

    for (int i=0; i<=N; ++i) {
        for (int j=0; j<(N+1)*(N+1); ++j) {
            Matrix[(Start_I + i*Inc)*(N+1)*(N+1) + j] = Matrix[(Start_I + Copy_I + i*Inc)*(N+1)*(N+1) + j];
        }
    }

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
void DG_BoundaryElement_2d::PeriodicBoundary(double *Matrix, initializer_list<int> I) {
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
void DG_BoundaryElement_2d::updateNeumann(string v, double *Matrix) {
    if (TopBoundary.count(v)) {
        if (TopBoundary[v] == "neumann") NeumannBoundary(Matrix, {N*(N+1), 1, -(N+1)});
    }
    if (BottomBoundary.count(v)) {
        if (BottomBoundary[v] == "neumann") NeumannBoundary(Matrix, {0, 1, (N+1)});
    }
    if (RightBoundary.count(v)) {
        if (RightBoundary[v] == "neumann") NeumannBoundary(Matrix, {N, N+1, -1});
    }
    if (LeftBoundary.count(v)) {
        if (LeftBoundary[v] == "neumann") NeumannBoundary(Matrix, {0, N+1, 1});
    }
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
void DG_BoundaryElement_2d::updateDirichlet(string v, double *Matrix) {
    if (TopBoundary.count(v)) {
        if (TopBoundary[v] == "dirichlet") DirichletBoundary(Matrix, {N*(N+1), 1});
    }
    if (BottomBoundary.count(v)) {
        if (BottomBoundary[v] == "dirichlet") DirichletBoundary(Matrix, {0, 1});
    }
    if (RightBoundary.count(v)) {
        if (RightBoundary[v] == "dirichlet") DirichletBoundary(Matrix, {N, N+1});
    }
    if (LeftBoundary.count(v)) {
        if (LeftBoundary[v] == "dirichlet") DirichletBoundary(Matrix, {0, N+1});
    }
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to get the variable `v` differentiated partially w.r.t. `x` and then store it in the
 * variable `vDash`. The function also takes `fluxType` as an input which would describe the numerical scheme that
 * should be used in order to obtain the derivative, modified for a boundary cell.
 *
 * @Param v         Variable which is to be differentiated.
 * @Param vDash     Variable in which the derivative is to be stored.
 * @Param fluxType  The type of flux that is to be used. eg "central"
 */
/* ----------------------------------------------------------------------------*/
void DG_BoundaryElement_2d::delByDelX(string v, string vDash, string conserVar, string fluxType, string fluxVariable = "") {
    double dy = (y_end - y_start);
    double dx = (x_end - x_start);

    double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
    double* RMatrix              =   new double[(N+1)*(N+1)*(N+1)*(N+1)]; // Creating a new temporary matrix, to store the modified RHS Matrices
    
    zeros(numericalFlux, (N+1)*(N+1));
    
    if(fluxType == "central") {    
                                              
        for(int i=0; i<=N; i++){
            numericalFlux[i*(N+1)+N]    = 0.5*( *boundaryRight[v][i]    + *neighboringRight[v][i] ) ;   
            numericalFlux[i*(N+1)]    = 0.5*( *boundaryLeft[v][i]     + *neighboringLeft[v][i] ) ;  
        }
        
    }

    else if(fluxType == "rusanov") {                                                    
        for(int i=0; i<=N; i++){
          // Normals nx, ny of the cell have been incorporated into the signs, need to set them separately !!
            numericalFlux[i*(N+1)+N] = 0.5*(*boundaryRight[v][i] + *neighboringRight[v][i] + MAX(fabs(*boundaryRight[fluxVariable][i]), fabs(*neighboringRight[fluxVariable][i]))*(*boundaryRight[conserVar][i] - *neighboringRight[conserVar][i])  ) ;   
            numericalFlux[i*(N+1)]   = 0.5*(*boundaryLeft[v][i]  + *neighboringLeft[v][i]  - MAX(fabs(*boundaryLeft[fluxVariable][i]), fabs(*neighboringLeft[fluxVariable][i]))*(*boundaryLeft[conserVar][i] - *neighboringLeft[conserVar][i])  ) ;   
        }

    }

    
    if (RightBoundary.count(v)) {
        if (RightBoundary[v] == "outflow") {
            for(int i=0; i<=N; i++){
          // Normals nx, ny of the cell have been incorporated into the signs, need to set them separately !!
            numericalFlux[i*(N+1)+N] = 0.5*(*boundaryRight[v][i] + *neighboringRight[v][i] + 2.0*fabs(*boundaryRight[fluxVariable][i])*(*boundaryRight[conserVar][i])  ) ;      
        }
        }
    }
    if (LeftBoundary.count(v)) {
        if (LeftBoundary[v] == "outflow") {
            for(int i=0; i<=N; i++){
          // Normals nx, ny of the cell have been incorporated into the signs, need to set them separately !!
            numericalFlux[i*(N+1)]   = 0.5*(*boundaryLeft[v][i]  + *neighboringLeft[v][i]  - 2.0*fabs(*boundaryLeft[fluxVariable][i])*(*boundaryLeft[conserVar][i])  ) ;   
        }
        }
    }

    // Derivative Matrix
    /// vDash = -0.5*dy*D*v
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, (N+1)*(N+1), (N+1)*(N+1), (N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix, (N+1)*(N+1), derivativeMatrix_x, (N+1)*(N+1), 0, RMatrix, (N+1)*(N+1));

    // Implementing Boundary Condition
    updateDirichlet(conserVar, RMatrix);
    updateNeumann(conserVar, RMatrix);
        
    cblas_dgemv(CblasRowMajor, CblasNoTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dy, RMatrix, (N+1)*(N+1), variable[v], 1, 0, variable[vDash], 1);

    /// Adding the numeical Flux terms as necessary.

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (N+1)*(N+1), (N+1)*(N+1), (N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix, (N+1)*(N+1), fluxMatrix_right, (N+1)*(N+1), 0, RMatrix, (N+1)*(N+1));

    // Implementing Boundary Condition
    updateDirichlet(conserVar, RMatrix);
    updateNeumann(conserVar, RMatrix);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dy, RMatrix, (N+1)*(N+1), numericalFlux, 1, 1, variable[vDash], 1);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (N+1)*(N+1), (N+1)*(N+1), (N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix, (N+1)*(N+1), fluxMatrix_left, (N+1)*(N+1), 0, RMatrix, (N+1)*(N+1));

    // Implementing Boundary Condition
    updateDirichlet(conserVar, RMatrix);
    updateNeumann(conserVar, RMatrix);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dy, RMatrix, (N+1)*(N+1), numericalFlux, 1, 1, variable[vDash], 1);

    delete[] numericalFlux;
    delete[] RMatrix;

    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to get the variable `v` differentiated partially w.r.t. `y` and then store it in the
 * variable `vDash`. The function also takes `fluxType` as an input which would describe the numerical scheme that
 * should be used in order to obtain the derivative, modified for boundary cells.
 *
 * @Param v         Variable which is to be differentiated.
 * @Param vDash     Variable in which the derivative is to be stored.
 * @Param fluxType  The type of flux that is to be used. eg "central"
 */
/* ----------------------------------------------------------------------------*/
void DG_BoundaryElement_2d::delByDelY(string v, string vDash, string conserVar, string fluxType, string fluxVariable = "") {
    double dy = (y_end - y_start);
    double dx = (x_end - x_start);

    double* numericalFlux        =   new double[(N+1)*(N+1)]; /// Creating a temporary new variable.
    double* RMatrix              =   new double[(N+1)*(N+1)*(N+1)*(N+1)]; // Creating a new temporary matrix, to store the modified RHS Matrices
    
    zeros(numericalFlux, (N+1)*(N+1));

    if(fluxType == "central") {                                          
        for(int i=0; i<=N; i++){
            numericalFlux[N*(N+1)+i]    = 0.5*( *boundaryTop[v][i]    + *neighboringTop[v][i] ) ;   
            numericalFlux[i]            = 0.5*( *boundaryBottom[v][i]     + *neighboringBottom[v][i] ) ;  
        }
        
    }
    
    else if(fluxType == "rusanov") {                                                     
        for(int i=0; i<=N; i++){
          // Normals nx, ny have been incorporated into the signs, need to set them separately !!
            numericalFlux[N*(N+1)+i]= 0.5*(*boundaryTop[v][i] + *neighboringTop[v][i] + MAX(fabs(*boundaryTop[fluxVariable][i]), fabs(*neighboringTop[fluxVariable][i]))*(*boundaryTop[conserVar][i] - *neighboringTop[conserVar][i]));
            numericalFlux[i]        = 0.5*(*boundaryBottom[v][i] + *neighboringBottom[v][i] - MAX(fabs(*boundaryBottom[fluxVariable][i]), fabs(*neighboringBottom[fluxVariable][i]))*(*boundaryBottom[conserVar][i] - *neighboringBottom[conserVar][i])); 
        }
      
    }

    if (TopBoundary.count(conserVar)) {
        if (TopBoundary[v] == "outflow") {
            for(int i=0; i<=N; i++){
            numericalFlux[N*(N+1)+i]    = 0.5*(*boundaryTop[v][i] + *neighboringTop[v][i]) + fabs(*boundaryTop[fluxVariable][i])*(*boundaryTop[conserVar][i]) ;    
            }
        }
    }
    if (BottomBoundary.count(v)) {
        if (BottomBoundary[v] == "outflow") {
            for(int i=0; i<=N; i++){ 
            numericalFlux[i]            = 0.5*(*boundaryBottom[v][i] + *neighboringBottom[v][i] - 2.0*fabs(*boundaryBottom[fluxVariable][i])*(*boundaryBottom[conserVar][i])) ;  
           }
        }
    }

    // Derivative Matrix
    /// vDash = -0.5*dy*D*v
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, (N+1)*(N+1), (N+1)*(N+1), (N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix, (N+1)*(N+1), derivativeMatrix_y, (N+1)*(N+1), 0, RMatrix, (N+1)*(N+1));

    // Implementing Boundary Condition
    updateDirichlet(conserVar, RMatrix);
    updateNeumann(conserVar, RMatrix);
        
    cblas_dgemv(CblasRowMajor, CblasNoTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dx, RMatrix, (N+1)*(N+1), variable[v], 1, 0, variable[vDash], 1);

    /// Adding the numeical Flux terms as necessary.

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (N+1)*(N+1), (N+1)*(N+1), (N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix, (N+1)*(N+1), fluxMatrix_top, (N+1)*(N+1), 0, RMatrix, (N+1)*(N+1));

    // Implementing Boundary Condition
    updateDirichlet(conserVar, RMatrix);
    updateNeumann(conserVar, RMatrix);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dx, RMatrix, (N+1)*(N+1), numericalFlux, 1, 1, variable[vDash], 1);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (N+1)*(N+1), (N+1)*(N+1), (N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix, (N+1)*(N+1), fluxMatrix_bottom, (N+1)*(N+1), 0, RMatrix, (N+1)*(N+1));

    // Implementing Boundary Condition
    updateDirichlet(conserVar, RMatrix);
    updateNeumann(conserVar, RMatrix);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dx, RMatrix, (N+1)*(N+1), numericalFlux, 1, 1, variable[vDash], 1);

    delete[] numericalFlux;
    delete[] RMatrix;

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
void DG_BoundaryElement_2d::limitMoments(string m, string modm, string cm, unsigned Index) {

 if (*variable[cm] && PositivityMarker) 
  { // Checking if cell marker is not equal to zero
    int count, Tempi, Tempj, i, j;
    count = N+1;
    double Temp1, Temp2, AlphaN;
    double epsilon = 1e-13;
    //AlphaN = sqrt((2.0*N -1.0)/(2.0*N +1)); // Similar to a diffusion coefficient
    //AlphaN = 0.5/sqrt(4.0*N*N -1.0);

    // Ensuring that Cell avergae remains the  same after limiting !!
    variable[modm][0] = variable[m][0];

    if ( N == 1) {
        epsilon = 0.0 ;
    }
    else {
        epsilon = 0.0;
    }

    for(i=Index; i > 0; i = i - (N+2)) {
     --count;
     //AlphaN = sqrt((2.0*(count)-1.0)/(2.0*(count)+1.0));
     AlphaN = 0.5*sqrt((4.0*(count)-1.0)/(2.0*(count)+1.0));
     //AlphaN = 0.5/sqrt(4.0*count*count -1.0);
     for(j=0; j < count; ++j) {
       Tempi = i-j;
       Tempj = i - j*(N+1);
       //Temp1 = BoundaryMinMod(m, Tempi, AlphaN, rightNeighbor, leftNeighbor, topNeighbor, bottomNeighbor);
       //Temp2 = BoundaryMinMod(m, Tempj, AlphaN, rightNeighbor, leftNeighbor, topNeighbor, bottomNeighbor);
       Temp1 = MinMod(variable[m][Tempi], AlphaN*(rightNeighbor->variable[m][Tempi-1] -variable[m][Tempi-1]), AlphaN*(variable[m][Tempi-1] -leftNeighbor->variable[m][Tempi-1]) , AlphaN*(topNeighbor->variable[m][Tempi-(N+1)] -variable[m][Tempi-(N+1)]), AlphaN*(variable[m][Tempi-(N+1)] -bottomNeighbor->variable[m][Tempi-(N+1)]));
       Temp2 = MinMod(variable[m][Tempj], AlphaN*(rightNeighbor->variable[m][Tempj-1] -variable[m][Tempj-1]), AlphaN*(variable[m][Tempj-1] -leftNeighbor->variable[m][Tempj-1]) , AlphaN*(topNeighbor->variable[m][Tempj-(N+1)] -variable[m][Tempj-(N+1)]), AlphaN*(variable[m][Tempj-(N+1)] -bottomNeighbor->variable[m][Tempj-(N+1)]));
       
       if ( abs(Temp1-variable[modm][Tempi]) > epsilon || abs(Temp2-variable[modm][Tempj]) > epsilon ) {
         variable[modm][Tempi] = Temp1;
         variable[modm][Tempj] = Temp2;
       }
       else //if(Temp1 != 0.0 && Temp2 != 0.0)
       {
         return ; // Need to exit both loops
       }
     }
     // Special Case for end values, when Tempi or Tempj access zero order polynomials !!
       Tempi = i-j;
       Tempj = i - j*(N+1);
       
       //Temp1 = BoundaryMinMod(m, Tempi, AlphaN, this, this, topNeighbor, bottomNeighbor);
       //Temp2 = BoundaryMinMod(m, Tempj, AlphaN, rightNeighbor, leftNeighbor, this, this);
       Temp1 = MinMod(variable[m][Tempi], AlphaN*(topNeighbor->variable[m][Tempi-(N+1)] -variable[m][Tempi-(N+1)]), AlphaN*(variable[m][Tempi-(N+1)] -bottomNeighbor->variable[m][Tempi-(N+1)]));
       Temp2 = MinMod(variable[m][Tempj], AlphaN*(rightNeighbor->variable[m][Tempj-1] -variable[m][Tempj-1]), AlphaN*(variable[m][Tempj-1] -leftNeighbor->variable[m][Tempj-1]));
      
       if ( abs(Temp1-variable[modm][Tempi]) > epsilon || abs(Temp2-variable[modm][Tempj]) > epsilon ) {
         variable[modm][Tempi] = Temp1;
         variable[modm][Tempj] = Temp2;
       }
       else //if(Temp1 != 0.0 && Temp2 != 0.0)
       {
         return ; // Need to exit both loops
       }
    }
 }

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
double DG_BoundaryElement_2d::BoundaryMinMod(string m, int Index, double Alpha, DG_Element_2d* R, DG_Element_2d* L, DG_Element_2d* T, DG_Element_2d* B) {
    vector<double> Elements;
    double epsilon = 1e-16;
    double M = 50.0;
    double min_dx = min(abs(X[0]-X[1]), abs(Y[0]-Y[1]));
    //Elements.push_back(M*min_dx*min_dx);

    Elements.push_back(variable[m][Index]);
    if (R != this)// && abs(R->variable[m][Index]) > epsilon ) 
    {
        Elements.push_back(Alpha*(R->variable[m][Index-1] -variable[m][Index-1]));
       // Elements.push_back((R->variable[m][Index]));
    } 
    if (L != this )//&& abs(L->variable[m][Index]) > epsilon )
     {
        Elements.push_back(Alpha*(variable[m][Index-1] - L->variable[m][Index-1]));
       // Elements.push_back((L->variable[m][Index]));
    } 
    if (T != this )//&& abs(T->variable[m][Index]) > epsilon )
     {
        Elements.push_back(Alpha*(T->variable[m][Index-(N+1)] -variable[m][Index-(N+1)]));
       // Elements.push_back((T->variable[m][Index]));
    } 
    if (B != this)// && abs(B->variable[m][Index]) > epsilon )
     {
        Elements.push_back(Alpha*(variable[m][Index-(N+1)] - B->variable[m][Index-(N+1)]));
        // Elements.push_back((B->variable[m][Index]));
    } 
    
    return MinMod(Elements);
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function updates Variables at the Boundary.
 *
 * @Param v This is the variable ,whose boundary values are to be fixed.
 */
/* ----------------------------------------------------------------------------*/
void DG_BoundaryElement_2d::updateBoundaryVariables(string v) {
    // Update only for Neumann/Outflow/Zero Gradient ,right now !!
    if (TopBoundary.count(v)) {
        if (TopBoundary[v] == "neumann") {
            for(int i=0; i<=N; ++i)
                variable[v][N*(N+1) + i] = variable[v][N*(N+1) + i -(N+1)];
        }
        /*else if (TopBoundary[v] == "dirichlet") {
            for(int i=0; i<=N; ++i)
                 *boundaryTop[v][i] = DirichletTop[v][i];
        } */ 
    }
    if (BottomBoundary.count(v)) {
        if (BottomBoundary[v] == "neumann") {
            for(int i=0; i<=N; ++i)
                variable[v][0 + i] = variable[v][0 + i +(N+1)];
        } 
       /* else if (BottomBoundary[v] == "dirichlet") {
            for(int i=0; i<=N; ++i)
                *boundaryBottom[v][i] = DirichletBottom[v][i];
        } */
    }
    if (RightBoundary.count(v)) {
        if (RightBoundary[v] == "neumann") {
            for(int i=0; i<=N; ++i)
                variable[v][N + i*(N+1)] = variable[v][N + i*(N+1) -1];
        } 
        /*else if (RightBoundary[v] == "dirichlet") {
            for(int i=0; i<=N; ++i)
                *boundaryRight[v][i] = DirichletRight[v][i];
        }*/
    }
    if (LeftBoundary.count(v)) {
        if (LeftBoundary[v] == "neumann"){
            for(int i=0; i<=N; ++i)
                variable[v][0 + i*(N+1)] = variable[v][0 + i*(N+1) +1];
        } 
        /*else if (LeftBoundary[v] == "dirichlet") {
            for(int i=0; i<=N; ++i)
                *boundaryLeft[v][i] = DirichletLeft[v][i];
        }*/
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
 * @Param cm This is the cell marker used to identified troubled cells.
*/
/* ----------------------------------------------------------------------------*/
void DG_BoundaryElement_2d::convertMomentToVariable(string m, string v, string cm) {
  /// Multiplying  VanderMand Matrix with the moments to obtained the nodal values of the variable.

 if (*variable[cm] && PositivityMarker)
  { 
      double *AuxVariable = new double[(N+1)*(N+1)];

      cblas_dgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 1.0, vanderMandMatrix,(N+1)*(N+1), variable[m],1,0,AuxVariable,1);

      if (TopBoundary.count(v)) {
        if (TopBoundary[v] == "neumann") {
            for(int i=0; i<=N; ++i)
                AuxVariable[N*(N+1) + i] = AuxVariable[N*(N+1) + i -(N+1)];
        }
        else if (TopBoundary[v] == "dirichlet") {
            for(int i=0; i<=N; ++i)
                 AuxVariable[N*(N+1) + i] = variable[v][N*(N+1) + i] ;
        }  
    }
    if (BottomBoundary.count(v)) {
        if (BottomBoundary[v] == "neumann") {
            for(int i=0; i<=N; ++i)
                AuxVariable[0 + i] = AuxVariable[0 + i +(N+1)];
        } 
        else if (BottomBoundary[v] == "dirichlet") {
            for(int i=0; i<=N; ++i)
                AuxVariable[0 + i] = variable[v][0 + i] ;
        } 
    }
    if (RightBoundary.count(v)) {
        if (RightBoundary[v] == "neumann") {
            for(int i=0; i<=N; ++i)
                AuxVariable[N + i*(N+1)] = AuxVariable[N + i*(N+1) -1];
        } 
        else if (RightBoundary[v] == "dirichlet") {
            for(int i=0; i<=N; ++i)
                AuxVariable[N + i*(N+1) ] = variable[v][N + i*(N+1)] ;
        }
    }
    if (LeftBoundary.count(v)) {
        if (LeftBoundary[v] == "neumann"){
            for(int i=0; i<=N; ++i)
                AuxVariable[0 + i*(N+1)] = AuxVariable[0 + i*(N+1) +1];
        } 
        else if (LeftBoundary[v] == "dirichlet") {
            for(int i=0; i<=N; ++i)
               AuxVariable[0 + i*(N+1)] = variable[v][0 + i*(N+1)] ;
        }
    }

    cblas_dcopy((N+1)*(N+1), AuxVariable, 1, variable[v], 1);

  /* for (int i=0; i< (N+1)*(N+1); ++i) {
      if (variable[v][i] < *variable["Min"]) {
          variable[v][i] = *variable["Min"] ;
      }
      else if (variable[v][i] > *variable["Max"]) {
          variable[v][i] = *variable["Max"] ;
      }
  }*/
    
    delete[] AuxVariable;
  

  }

  return ;
}



