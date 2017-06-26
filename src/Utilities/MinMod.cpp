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
  if ( Signum(a) == Signum(b) && Signum(b) == Signum(c) ) return (MIN(fabs(a), fabs(b), fabs(c))) * Signum(a);
  
  return 0.0 ;
}

double MinMod(double a, double b, double c, double d, double e) {
  if ( Signum(a) == Signum(b) && Signum(b) == Signum(c) && Signum(c) == Signum(d) && Signum(d) == Signum(e) ) 
    return (MIN(fabs(a), fabs(b), fabs(c), fabs(d), fabs(e))) * Signum(a);

  return 0.0 ;
}
