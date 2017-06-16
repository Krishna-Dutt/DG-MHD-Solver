#include "../../includes/Utilities/VanderMandLegendre.h"
#include "../../includes/Utilities/HeaderFiles.h"
#include "../../includes/Utilities/LegendrePolynomial.h"
#include "../../includes/Utilities/PolyEval.h"
#include "../../includes/Utilities/LobattoNodes.h"

using namespace std;

void twoDVanderMandLegendre(double *VanderMandMatrix, unsigned N) {
  double *Nodes = new double [N+1];
  double **LegendrePoly;

  LegendrePoly = new double*[N+1];
  lobattoNodes(Nodes, N+1);

  unsigned i1,j1,i2,j2;
  double k;

  for( i=0; i<=N; ++i) {
    LegendrePoly[i] = new double [i+1];
    legendrePolynomial(LegendrePoly[i], i);
  }

   function<double(double)> eval;
   for(i1=0;i1<=N;i1++)
    {
        for(j1=0;j1<=N;j1++)
        {
          for(i2=0; i2<=N; ++i2)
            for(j2=0 ; j2<=N; ++j2) {
              k = sqrt((2.0*i2+1.0) * (2.0*j2+1.0))*0.5;
              VanderMandMatrix[i*(N+1)+j] = k * polyEval(LegendrePoly[i],i,Nodes[i]) * polyEval(LegendrePoly[j],j,Nodes[j]);
        }
    }

   delete[] Nodes;
   for(i=0; i<=N; ++i) {
     delete[] LegendrePoly[i];
   }
   delete[] LegendrePoly;

   return ;
}


 


