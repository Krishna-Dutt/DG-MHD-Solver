#include "../../includes/Utilities/VanderMandLegendre.h"
#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/LegendrePolynomial.h"
#include "../../includes/Utilities/PolyEval.h"
#include "../../includes/Utilities/LobattoNodes.h"

using namespace std;

void VanderMandLegendre(double *VanderMandMatrix, unsigned N) {
  double *Nodes = new double [N+1];
  double **LegendrePoly;

  LegendrePoly = new double*[N+1];
  lobattoNodes(Nodes, N+1);

  unsigned i,j;

  for( i=0; i<=N; ++i) {
    LegendrePoly[i] = new double [i+1];
    legendrePolynomial(LegendrePoly[i], i);
  }

   for(i=0;i<=N;i++)
    {
      for(j=0;j<=N;j++) {
              VanderMandMatrix[i*(N+1)+j] = polyEval(LegendrePoly[j],j,Nodes[i]);
        }
    }

   delete[] Nodes;
   for(i=0; i<=N; ++i) {
     delete[] LegendrePoly[i];
   }
   delete[] LegendrePoly;

   return ;
}


void twoDVanderMandLegendre(double *VanderMandMatrix, unsigned N) {
  double VM[(N+1)][(N+1)];

  int i1,i2,j1,j2;
  double k;

  VanderMandLegendre(*VM, N);

  for(i1=0; i1<=N; ++i1)
    for(j1=0; j1<=N; ++j1) {
      for(i2=0; i2<=N; ++i2)
        for(j2=0; j2<=N; ++j2) {
          k = sqrt((2.0*i2+1.0) * (2.0*j2+1.0)) * 0.5;
          VanderMandMatrix[(i1*(N+1)+j1)*(N+1)*(N+1)+i2*(N+1)+j2] = k * VM[i1][i2] * VM[j1][j2];
        }
    }
  return ;
}

 
 


