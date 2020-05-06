#include <cmath>
#include "PoissonFunction.hpp"


// Constructor
PoissonFunction::PoissonFunction()
{}

// Destructor
PoissonFunction::~PoissonFunction()
{}

// Overriding pure virtual methods
double PoissonFunction::evaluateF(double x, double y)
{
  return (pow(M_PI,2.0)/4.0) * cos(0.5*M_PI*x) * (1.0 - 17.0*cos(2.0*M_PI*y));
}

double PoissonFunction::evaluateBoundary(double x, double y)
{
  return cos(0.5*M_PI*x) * (1 - cos(2.0*M_PI*y));
}
