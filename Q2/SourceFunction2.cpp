#include <cmath>
#include "SourceFunction2.hpp"


// Constructor
SourceFunction2::SourceFunction2()
{}

// Destructor
SourceFunction2::~SourceFunction2()
{}

// Overriding pure virtual methods
double SourceFunction2::evaluateF(double x, double y)
/* Evaluates f(x,y) */
{
  return (pow(M_PI,2.0)/4.0) * cos(0.5*M_PI*x) * (1.0 - 17.0*cos(2.0*M_PI*y));
}


double SourceFunction2::evaluateBoundary(double x, double y)
/* Evaluates g(x,y), i.e. boundary function */
{
  return cos(0.5*M_PI*x) * (1.0 - cos(2.0*M_PI*y));
}


double SourceFunction2::evaluateUexact(double x, double y)
/* Evaluates u_exact(x,y) */
{
  return cos(0.5*M_PI*x) * (1.0 - cos(2.0*M_PI*y));
}
