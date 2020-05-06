#ifndef POISSONFUNCTION
#define POISSONFUNCTION

#include "Abstract2DFunction.hpp"

class PoissonFunction: public Abstract2DFunction
{
  public:

    // Constructor
    PoissonFunction();

    // Destructor
    ~PoissonFunction();

    // Overriding pure virtual methods
    double evaluateF(double x, double y);
    double evaluateBoundary(double x, double y);

};

#endif
