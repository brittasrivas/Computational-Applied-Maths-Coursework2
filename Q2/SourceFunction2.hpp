#ifndef SOURCEFUNCTION2
#define SOURCEFUNCTION2

#include "Abstract2DFunction.hpp"

class SourceFunction2: public Abstract2DFunction
{
  public:

    // Constructor
    SourceFunction2();

    // Destructor
    ~SourceFunction2();

    // Overriding pure virtual methods
    double evaluateF(double x, double y);
    double evaluateBoundary(double x, double y);

};

#endif
