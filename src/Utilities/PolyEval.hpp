#include <lapacke.h>

#ifndef PolyEval_HPP
#define PolyEval_HPP

float polyEval(float *GivenPoly, unsigned deg, float x)
{
	float Result = GivenPoly[deg];
	for(int i = deg-1;i>=0;i--)
		Result = Result*x+GivenPoly[i];
	return Result ;
}

#endif