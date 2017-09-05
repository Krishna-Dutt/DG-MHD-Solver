#include "../../includes/Utilities/MinMod.h"
#include "../../includes/Utilities/HeaderFiles.h"

double MIN(double x, double y, double z) {
  return min(min(x,y), min(y,z)) ;
}

double MIN(double x, double y, double z, double a, double b) {
  return min(MIN(x,y,z), min(a,b)) ;
}

double Signum(double x) {
  if ( x > 0.0) return 1.0 ;
  if ( x < 0.0) return -1.0 ;
  return 0.0 ;
}

double MinMod(double a, double b, double c) {
 //if ( Signum(a) == Signum(b) && Signum(b) == Signum(c) ) 
 // if ( Signum(b) == Signum(c) ) 
  return (MIN(fabs(a), fabs(b), fabs(c))) * Signum(a);
  
  return 0.0 ;
}

double MinMod(double a, double b, double c, double d, double e) {
  //if ( Signum(a) == Signum(b) && Signum(b) == Signum(c) && Signum(c) == Signum(d) && Signum(d) == Signum(e) ) 
  //if (  Signum(b) == Signum(c) && Signum(d) == Signum(e) ) 
    return (MIN(fabs(a), fabs(b), fabs(c), fabs(d), fabs(e))) * Signum(a);

  return 0.0 ;
}

// Implementing Modified version of Minmod that works for N=1 Limiting !!
double MinMod(vector<double> A) {
  double Sign = Signum(A[0]);

  double Min = abs(A[0]);
  for(int i=1; i<A.size(); ++i) {
    Min = min(abs(A[i]), Min);
    if (Signum(A[i]) != Sign) return 0.0 ;
  }

  return Min*Signum(A[0]);
  /*// TVB Minmod 
  double Sign =  Signum(A[1]);
  if (abs(A[1]) <= abs(A[0])) return A[1];

  double Min = abs(A[1]);
  for(int i=1; i<A.size(); ++i) {
    Min = min(abs(A[i]), Min);
    //if (Signum(A[i]) != Sign) return 0.0 ;
  }

  return Min*Signum(A[1]);*/
  /*double Min = abs(A[0]);
  int Sign = Signum(A[0]);
  int size = A.size();

  if (size == 3) {
     //if (Signum(A[1]) == Signum(A[2]))
      return  Sign*MIN(A[0],min(2.0*A[1],A[2]),min(A[1], 2.0*A[2]));
     return 0.0 ;
  }
  else if (size == 5) {
     //if (Signum(A[1]) == Signum(A[2]) && Signum(A[3]) == Signum(A[4])) 
     {
       return  Sign*MIN(A[0],min(2.0*A[1],A[2]),min(A[1], 2.0*A[2]),min(2.0*A[3],A[4]),min(A[3], 2.0*A[4]));
     }
     return 0.0 ;

  }*/

}