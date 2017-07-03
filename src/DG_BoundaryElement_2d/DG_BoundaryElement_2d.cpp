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
            break;
        case 'r' : // `r` or `R` for right
        case 'R' :
            RightBoundary[v] = type;
            break;
        case 'b' : // `b` or `B` for bottom
        case 'B' :
            BottomBoundary[v] = type;
            break;
        case 'l' : // `l` or `L` for left
        case 'L' :
            LeftBoundary[v] = type;
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